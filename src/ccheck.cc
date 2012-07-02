#include <map>
#include <string>
#include <utility>

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>

extern "C" {
#include "map_align.h"
#include "mia.h"
#include "myers_align.h"
}

/*
 * Contamination Checker.  Outline:
 *
 * - read the human reference; this will contain ambiguity codes
 * - read maln file, including assembly and assembled reads
 * - align reference-consensus and assembly globally
 *   This uses Myers' O(nd) aligner, for it grasps ambiguity codes and
 *   runs fast enough for long, but similar sequences.
 * - find "diagnostic positions", positions where ass and ref differ
 * - for every "end" fragment: store it  and later join with its other
 *   half
 * - for every "full" fragment: if it crosses at least one diagnostic
 *   position, cut out that range from ref and align to it globally
 *   using the mia aligner
 * - for every position where the bases agree, classify it, then
 *   classify the fragment (conflicting, uninformative, contaminant,
 *   endogenous)
 * - produce a summary
 *
 * Notable features:
 * - uses substitution matrix from maln file
 * - operates sensibly on aDNA
 * - has sensible commandline and doesn't make too much noise in operation
 * - optionally considers only certain diagnostic positions
 *   (tranversions only and/or some region only)
 *
 * To be done:
 * - read multiple maln files if desired and sum up statistics
 *   [not needed right now and postponed]
 * - filter out some diagnostic positions depending on shape of sequence id
 *   [not needed right now and postponed]
 * - filter out some differences as not being diagnostic (esp. near
 *   homopolymers; see existing contamination checker)
 *   [low priority for now and postponed]
 *
 * Nice to have:
 * - check if small/capital letters and all ambiguity codes work
 *   everywhere sensibly, then make use of them
 * - somehow account for sequencing error, especially on Solexa, when
 *   calculating the likely contamination level
 * - remove linear scans where binary search methods would work
 */

void print_aln( const char* aln1, const char* aln2 )
{
	int p ;
	const char* a ;
	while( *aln1 && *aln2 ) {
		for( p = 0, a = aln1 ; *a && p != 72 ; ++p ) putc( *a++, stderr ) ;
		putc( '\n', stderr ) ;

		for( p = 0, a = aln2 ; *a && p != 72 ; ++p ) putc( *a++, stderr ) ;
		putc( '\n', stderr ) ;

		for( p = 0 ; *aln1 && *aln2 && p != 72 ; ++p )
			putc( *aln1++ == *aln2++ ? '*' : ' ', stderr ) ;
		putc( '\n', stderr ) ;
		putc( '\n', stderr ) ;
	}
}

// List of diagnostic positions: Coordinates are relative to assembly
// (we want to quickly know whether a fragment overlaps a DP).  We'll
// store the reference bases along with it.
typedef std::map< int, std::pair< char, char > > dp_list ;

// Everything that differs counts as diagnostic, unless it's a gap.  In
// principle, Ns could be diagnostic, too, even though in only one
// direction.  In practice, however, it turned out that Ns produce noise
// and little in the way of useable results.  So Ns don't count as
// diagnostic for now.
bool is_diagnostic( const char* aln1, const char* aln2 )
{
	return *aln1 != *aln2
				&& *aln1 != 'N' && *aln2 != 'N' 
				&& *aln1 != '-' && *aln2 != '-' ;
}

bool is_transversion( char a, char b )
{
	char u = a & ~32 ;
	char v = b & ~32 ;
	switch( u )
	{
		case 'A': return v != 'G' ;
		case 'C': return v != 'T' ;
		case 'G': return v != 'A' ;
		case 'T':
		case 'U': return v != 'C' ;
		default: return false ;
	}
}


dp_list mk_dp_list( const char* aln1, const char* aln2, bool transversions, int span_from, int span_to, int& total_d )
{
	dp_list l ;
    int index = 0 ;
    total_d = 0 ;
    while( index != span_from && *aln1 && *aln2 )
    {
		if( *aln2 != '-' ) ++index ;
		++aln1 ;
		++aln2 ;
    }
	while( index != span_to && *aln1 && *aln2 )
	{
        if( *aln1 != *aln2 && *aln1 != '-' && *aln2 != '-' ) ++total_d ;
		if( is_diagnostic( aln1, aln2 ) && ( !transversions || is_transversion( *aln1, *aln2 )))
			l[index] = std::make_pair( *aln1, *aln2 ) ;
		if( *aln2 != '-' ) ++index ;
		++aln1 ;
		++aln2 ;
	}
	return l ;
}

std::pair< dp_list::const_iterator, dp_list::const_iterator >
overlapped_diagnostic_positions( const dp_list& l, const AlnSeqP s )
{
	dp_list::const_iterator left  = l.lower_bound( s->start ) ;
	dp_list::const_iterator right = l.lower_bound( s->end + 1 ) ;
	return std::make_pair( left, right ) ;
}

// XXX: linear scan --> O(n)
// This could be faster (O(log n)) if precompiled into some sort of index.
std::string lift_over( const char* aln1, const char* aln2, int s, int e ) 
{
	std::string r ;
	int p ;
	for( p = 0 ; p < e && *aln1 && *aln2 ; ++aln1, ++aln2 )
	{
		if( *aln1 != '-' && p >= s ) r.push_back( *aln1 ) ;
		if( *aln2 != '-' ) ++p ;
	}
	return r ;
}
 
bool consistent( bool adna, char x, char y )
{
	char x_ = x == 'G' ? 'R' : x == 'C' ? 'Y' : x ;
	return x == '-' || y == '-' || (char_to_bitmap( adna ? x_ : x ) & char_to_bitmap(y)) != 0 ;
}

bool check_perfect_match_1( bool adna, const std::string& ass, const char* seq )
{
    std::string::const_iterator c = ass.begin(), ce = ass.end() ; 
    while( c != ce && *seq )
    {
        if( *c == 'N' || *seq == 'N' ) return false ;
        if( *c != *seq &&
            (!adna || *c != 'C' || *seq != 'T') &&
            (!adna || *c != 'G' || *seq != 'A') ) return false ;
        ++c ;
        ++seq ;
    }
    return( c == ce && !*seq ) ;
}

bool check_perfect_match( bool adna, const std::string& ass, const char* seq )
{
    bool b = check_perfect_match_1( adna, ass, seq ) ;
    fprintf( stderr, "%s\n%s -- %s\n", ass.c_str(), seq, b ? "match" : "mismatch" ) ;
    return b ;
}

enum whatsit { unknown, clean, dirt, conflict, nonsense, maxwhatsits } ;

const char *label[] = { "unclassified", "clean       ", "polluting   ", "conflicting ", "nonsensical " } ;

whatsit merge_whatsit( whatsit a, whatsit b )
{
	if( a == b ) return a ;
	if( a == unknown ) return b ;
	if( b == unknown ) return a ;
	if( a == nonsense || b == nonsense ) return nonsense ;
	return conflict ;
}

bool sanity_check_sequence( const char* s )
{
	for( ; *s ; ++s )
		if( !strchr( "ACGTUN", *s ) )
			return false ;
	return true ;
}

std::string find_maln( std::string fn )
{
    if( fn.length() <= 2 || fn.substr( fn.length()-2 ) != ".1" ) return fn ;
    size_t p = fn.rfind( '/' ) ;
    std::string dir = p == std::string::npos ? std::string(".") : fn.substr( 0, p ) ;
    std::string base = p == std::string::npos ? fn.substr( 0, fn.length()-1 ) 
                                              : fn.substr( p+1, fn.length()-p-2 ) ;
    int num = 1 ;
    DIR* d = opendir( dir.c_str() ) ;
    while( struct dirent* de = readdir( d ) )
    {
        if( strlen(de->d_name) > base.length() 
                && base == std::string( de->d_name, base.length() ) )
        {
            char *n1 = de->d_name + base.length(), *n2 = n1 ;
            while( *n1 && isdigit(*n1) ) ++n1 ;
            if( !*n1 ) {
                int n = atoi( n2 ) ;
                if( n > num ) {
                    num = n ;
                    fn = p == std::string::npos ? de->d_name 
                                                : dir + "/" + de->d_name ;
                }
            }
        }
    }
    closedir( d ) ;
    return fn ;
}

struct option longopts[] = {
	{ "reference", required_argument, 0, 'r' },
	{ "ancient", no_argument, 0, 'a' },
	{ "verbose", no_argument, 0, 'v' },
	{ "help",    no_argument, 0, 'h' },
	{ "transversions", no_argument, 0, 't' },
	{ "span", required_argument, 0, 's' },
	{ "maxd", required_argument, 0, 'd' },
    { "table", no_argument, 0, 'T' },
	{ 0,0,0,0 }
} ;

void usage( const char* pname )
{
	fputs( "Usage: ", stdout ) ;
	fputs( pname, stdout ) ;
	fputs( " [-r <ref.fa>] [-a] [-t] [-s M-N] [-v] <aln.maln> \n\n"
		"Reads a maln file and tries to quantify contained contamination.\n"
		"Options:\n"
		"  -r, --reference FILE     FASTA file with the likely contaminant (default: builtin mt311)\n"
		"  -a, --ancient            Treat DNA as ancient (i.e. likely deaminated)\n"
		"  -t, --transversions      Treat only transversions as diagnostic\n"
		"  -s, --span M-N           Look only at range from M to N\n"
		"  -n, --numpos N           Require N diagnostic sites in a single read (default: 1)\n"
        "  -M, --max-maln           Find the highest numbered maln file in a series\n"
        "  -T, --table              Output as tables (easier for scripts, herder on the eyes)\n"
		"  -v, --verbose            Increase verbosity level (can be repeated)\n"
		"  -h, --help               Print this help message\n\n", stdout ) ;
}

extern       char mt311_sequence[] ;
extern const int  mt311_sequence_size ;

struct refseq hum_ref = {
	"mt311", "consensus of 311 human mitochondria",
	mt311_sequence, NULL,
	mt311_sequence_size-1, mt311_sequence_size,
	NULL, 1, 0 } ;

int main( int argc, char * const argv[] )
{
	bool adna = false ;
	bool transversions = false ;
    bool be_clever = false ;
    bool mktable = false ;
	int min_diag_posns = 1 ;
	int verbose = 0 ;
	int maxd = 0 ;
	int span_from = 0, span_to = INT_MAX ;

	if( argc == 0 ) { usage( argv[0] ) ; return 0 ; }

	int opt ;
	do {
		opt = getopt_long( argc, argv, "r:avhts:d:n:MT", longopts, 0 ) ;
		switch( opt ) 
		{
			case 'r': 
				read_fasta_ref( &hum_ref, optarg ) ;
				break ;
			case 'a':
				adna = true ;
				break ;
			case 'v':
				++verbose ;
				break ;
			case ':':
				fputs( "missing option argument\n", stderr ) ;
				break ;
			case '?':
				fputs( "unknown option\n", stderr ) ;
				break ;
			case 'h':
				usage( argv[0] ) ;
				return 1 ;
			case 't':
				transversions = true ;
				break ;
			case 's':
				sscanf( optarg, "%u-%u", &span_from, &span_to ) ;
				if( span_from ) span_from-- ;
				break ;
			case 'n':
				min_diag_posns = atoi( optarg ) ;
				break ;
			case 'd':
				maxd = atoi( optarg ) ;
				break ;
            case 'M':
                be_clever = true ;
                break ;
            case 'T':
                mktable = true ;
                break ;
		}
	} while( opt != -1 ) ;

	if( optind == argc ) { usage( argv[0] ) ; return 1 ; }
	if( !hum_ref.rcseq ) make_reverse_complement( &hum_ref ) ;
	bool hum_ref_ok = sanity_check_sequence( hum_ref.seq ) ;
	if( !hum_ref_ok ) fputs( "FUBAR'ed FastA file: contaminant sequence contains gap symbols.\n", stderr ) ;

    if( mktable ) {
        fputs( "#Filename\tAln.dist\t#diff\t#diag\t#tv\t", stdout ) ;
        for( whatsit klass = unknown ; klass != maxwhatsits ; klass = (whatsit)( (int)klass +1 ) )
        {
            fputs( label[klass], stdout ) ;
            putchar( '\t' ) ;
        }
        puts( "LB\tML\tUB\t#match\t#mism\tLB'\tML'\tUB'" ) ;
    }

    for( ; optind != argc ; ++optind )
    {
        int summary[ maxwhatsits ] = {0} ;
        int matched = 0, mismatched = 0 ;

        std::string infile( be_clever ? find_maln( argv[optind] ) : argv[optind] ) ;
        if( mktable ) {
            fputs( infile.c_str(), stdout ) ;
            putchar( '\t' ) ;
        }
        else {
            puts( infile.c_str() ) ;
            putchar( '\n' ) ;
        }
        MapAlignmentP maln = read_ma( infile.c_str() ) ;
        PSSMP submat = maln->fpsm ;

        bool maln_ref_ok = sanity_check_sequence( maln->ref->seq ) ;
        if( !maln_ref_ok ) fputs( "FUBAR'ed maln file: consensus sequence contains gap symbols.\n", stderr ) ;
        if( !hum_ref_ok || !maln_ref_ok ) {
            fputs( "Problem might exist between keyboard and chair.  I give up.\n", stderr ) ;
            return 1 ;
        }

        if( !maxd ) maxd = max( strlen(hum_ref.seq), strlen(maln->ref->seq) ) / 10 ;
        char *aln_con = (char*)malloc( strlen(hum_ref.seq) + maxd + 2 ) ;
        char *aln_ass = (char*)malloc( strlen(maln->ref->seq) + maxd + 2 ) ;
        unsigned d = myers_diff( hum_ref.seq, myers_align_globally, maln->ref->seq, maxd, aln_con, aln_ass ) ;

        if( d == UINT_MAX ) {
            fprintf( stderr, " *** Could not align references with up to %d mismatches.\n"
                             " *** This is usually a sign of trouble, but\n"
                             " *** IF AND ONLY IF YOU KNOW WHAT YOU ARE DOING, you can\n"
                             " *** try the -d N option with N > %d.\n", maxd, maxd ) ;
            return 1 ;
        }
        if( mktable ) printf( "%d\t", d ) ;
        else printf( "  %d alignment distance between reference and assembly.\n", d ) ;

        if( verbose >= 6 ) print_aln( aln_con, aln_ass ) ;

        int total_dist ;
        dp_list l = mk_dp_list( aln_con, aln_ass, transversions, span_from, span_to, total_dist ) ;
        if( mktable ) printf( "%d\t", total_dist ) ;
        else printf( "  %d total differences between reference and assembly.\n", total_dist ) ;

        {
            int t = 0 ;
            for( dp_list::const_iterator i = l.begin() ; i != l.end() ; ++i )
                if( is_transversion( i->second.first, i->second.second ) ) ++t ;
            if( mktable ) printf( "%d\t%d\t", (int)l.size(), t ) ; 
            else {
                printf( "  %d diagnostic positions", (int)l.size() ) ;
                if( span_from != 0 || span_to != INT_MAX )
                    printf( " in range [%d,%d)", span_from, span_to ) ;
                printf( ", %d of which are transversions.\n\n", t ) ;
            }
        }
        if( verbose >= 3 ) 
        {
            dp_list::const_iterator i = l.begin() ;
            if( i != l.end() ) { fprintf( stderr, "<%d:%c,%c>", i->first, i->second.first, i->second.second ) ; ++i ; }
            for( ; i != l.end() ; ++i ) fprintf( stderr, ", <%d:%c,%c>", i->first, i->second.first, i->second.second ) ;
            putc( '\n', stderr ) ;
        }

        typedef std::map< std::string, std::pair< std::pair< whatsit, int >, bool > > Bfrags ;
        Bfrags bfrags ;
        const AlnSeqP *s ;

        for( s = maln->AlnSeqArray ; s != maln->AlnSeqArray + maln->num_aln_seqs ; ++s )
        {
            /* Fixup stupid naming in maln files.  Must everything be
             * encoded in the name, dammit?! */
            {
                char *p = (*s)->id, *q = (*s)->id ;
                while( *q ) ++q ;
                if( q-p > 3 && (q[-1] == 'b' || q[-1] == 'f') && q[-2] == '_' )
                {
                    if( q[-3] == ',' ) q[-3] = 0 ; else q[-2] = 0 ;
                }
            }

            whatsit klass = unknown ;
            int votes = 0 ;

            std::string the_ass( maln->ref->seq + (*s)->start, (*s)->end - (*s)->start + 1 ) ;
            bool is_perfect_match = verbose >= 5 ? check_perfect_match( adna, the_ass, (*s)->seq )
                                                 : check_perfect_match_1( adna, the_ass, (*s)->seq ) ;
            if( verbose >= 3 ) {
                fputs( (*s)->id, stderr ) ;
                putc( '/', stderr ) ;
                putc( (*s)->segment, stderr ) ;
                if( is_perfect_match ) fputs( ": perfect match\n", stderr ) ;
                else fputs( ": contains mismatches\n", stderr ) ;
            }


            std::pair< dp_list::const_iterator, dp_list::const_iterator > p =
                overlapped_diagnostic_positions( l, *s ) ;
            if( std::distance( p.first, p.second ) < min_diag_posns )
            {
                if( verbose >= 3 ) {
                    fputs( (*s)->id, stderr ) ;
                    putc( '/', stderr ) ;
                    putc( (*s)->segment, stderr ) ;
                    fputs( ": no diagnostic positions\n", stderr ) ;
                }
            }
            else
            {
                if( verbose >= 3 )
                {
                    fprintf( stderr, "%s/%c: %d diagnostic positions", (*s)->id, (*s)->segment, (int)std::distance( p.first, p.second ) ) ;
                    if( verbose >= 4 ) 
                    {
                        putc( ':', stderr ) ; putc( ' ', stderr ) ;
                        dp_list::const_iterator i = p.first ;
                        if( i != p.second ) { fprintf( stderr, "<%d:%c,%c>", i->first, i->second.first, i->second.second ) ; ++i ; }
                        for( ; i != p.second ; ++i ) fprintf( stderr, ", <%d:%c,%c>", i->first, i->second.first, i->second.second ) ;
                    }
                    fprintf( stderr, "\nrange:  %d..%d\n", (*s)->start, (*s)->end ) ;
                }

                std::string the_read ;
                for( char *nt = (*s)->seq, **ins = (*s)->ins ; *nt ; ++nt, ++ins )
                {
                    if( *nt != '-' ) the_read.push_back( *nt ) ;
                    if( *ins ) the_read.append( *ins ) ;
                }
                std::string lifted = lift_over( aln_con, aln_ass, (*s)->start, (*s)->end + 1 ) ;

                if( verbose >= 5 )
                {
                    fprintf( stderr, "raw read: %s\nlifted:   %s\nassembly: %s\n\n"
                                     "aln.read: %s\naln.assm: %s\nmatches:  ",
                                     the_read.c_str(), lifted.c_str(), the_ass.c_str(), 
                                     (*s)->seq, the_ass.c_str() ) ;
                    std::string::const_iterator b = the_ass.begin(), e = the_ass.end() ;
                    const char* pc = (*s)->seq ;
                    while( b != e && *pc ) putc( *b++ == *pc++ ? '*' : ' ', stderr ) ;
                    putc( '\n', stderr ) ;
                }

                int size = std::max( lifted.size(), the_read.size() ) ;

                AlignmentP frag_aln = init_alignment( size, size, 0, 0 ) ;

                frag_aln->seq1 = lifted.c_str() ;
                frag_aln->seq2 = the_read.c_str() ;
                frag_aln->len1 = size ;
                frag_aln->len2 = size ;
                frag_aln->sg5 = 1 ;
                frag_aln->sg3 = 1 ;
                frag_aln->submat = submat ;
                pop_s1c_in_a( frag_aln ) ;
                pop_s2c_in_a( frag_aln ) ;
                dyn_prog( frag_aln ) ;

                pw_aln_frag pwaln ;
                max_sg_score( frag_aln ) ;			// ARGH!  This has a vital side-effect!!!
                find_align_begin( frag_aln ) ;  	//        And so has this...
                populate_pwaln_to_begin( frag_aln, &pwaln ) ;
                pwaln.start = frag_aln->abc;

                if( verbose >= 5 )
                {
                    fprintf( stderr, "\naln.read: %s\naln.ref:  %s\nmatches:  ", pwaln.frag_seq, pwaln.ref_seq ) ;
                    const char *pc = pwaln.frag_seq, *pd = pwaln.ref_seq ;
                    while( *pc && *pd ) putc( *pd++ == *pc++ ? '*' : ' ', stderr ) ;
                    putc( '\n', stderr ) ;
                    putc( '\n', stderr ) ;
                }

                free_alignment( frag_aln ) ;

                char *paln1 = aln_con, *paln2 = aln_ass ;
                int ass_pos = 0 ;
                while( ass_pos != (*s)->start && *paln1 && *paln2 ) 
                {
                    if( *paln2 != '-' ) ass_pos++ ;
                    ++paln1 ;
                    ++paln2 ;
                }

                std::string in_ref = lifted.substr( 0, pwaln.start ) ;
                in_ref.append( pwaln.ref_seq ) ;

                char *in_frag_v_ref = pwaln.frag_seq ;
                char *in_ass = maln->ref->seq + (*s)->start ;
                char *in_frag_v_ass = (*s)->seq ;

                if( verbose ) {
                    if(*paln1!=in_ref[0]||*paln1=='-') fprintf( stderr, "huh? (R+%d) %.10s %.10s\n", pwaln.start, paln1, in_ref.c_str() ) ;
                    if(*paln2!=in_ass[0]&&*paln2!='-') fprintf( stderr, "huh? (A+%d) %.10s %.10s\n", pwaln.start, paln2, in_ass ) ;
                }

                while( ass_pos != (*s)->end +1 && *paln1 && *paln2 && !in_ref.empty() && *in_ass && *in_frag_v_ass && *in_frag_v_ref )
                {
                    if( is_diagnostic( paln1, paln2 ) ) {
                        if( verbose >= 4 )
                            fprintf( stderr, "diagnostic pos.: %d %c/%c %c/%c",
                                     ass_pos, in_ref[0], *in_frag_v_ref, *in_ass, *in_frag_v_ass ) ;
                        if( *in_frag_v_ref != *in_frag_v_ass ) 
                        {
                            if( verbose >= 4 ) fputs( "in disagreement.\n", stderr ) ;
                        }
                        else
                        {
                            bool maybe_clean = consistent( adna, *in_ass, *in_frag_v_ass ) ;
                            bool maybe_dirt =  consistent( adna, in_ref[0], *in_frag_v_ref ) ;

                            if( verbose >= 4 )
                            {
                                fputs( maybe_dirt  ? "" : "in", stderr ) ;
                                fputs( " consistent/", stderr ) ;
                                fputs( maybe_clean ? "" : "in", stderr ) ;
                                fputs( " consistent\n", stderr ) ; 
                            }

                            if( maybe_clean && !maybe_dirt && klass == unknown ) klass = clean ;
                            if( maybe_clean && !maybe_dirt && klass == dirt    ) klass = conflict ;
                            if( !maybe_clean && maybe_dirt && klass == unknown ) klass = dirt ;
                            if( !maybe_clean && maybe_dirt && klass == clean   ) klass = conflict ;
                            if( !maybe_clean && !maybe_dirt )                    klass = nonsense ;
                            if( maybe_dirt != maybe_clean ) votes++ ;
                        }
                    }

                    if( *paln1 != '-' ) {
                        do {
                            in_ref=in_ref.substr(1) ;
                            in_frag_v_ref++ ;
                        } while( in_ref[0] == '-' ) ;
                    }
                    if( *paln2 != '-' ) {
                        ass_pos++ ;
                        do {
                            in_ass++ ;
                            in_frag_v_ass++ ;
                        } while( *in_ass == '-' ) ;
                    }
                    ++paln1 ;
                    ++paln2 ;
                }
            }

            Bfrags::const_iterator i = bfrags.find( (*s)->id ) ;


            switch( (*s)->segment )
            {
                case 'b':
                    bfrags[ (*s)->id ] = std::make_pair( std::make_pair( klass, votes ), is_perfect_match ) ;
                    if( verbose >= 3 ) putc( '\n', stderr ) ;
                    break ;

                case 'f':
                    if( i == bfrags.end() ) 
                    {
                        fputs( (*s)->id, stderr ) ;
                        fputs( "/f is missing its back.\n", stderr ) ;
                    }
                    else
                    {
                        votes += i->second.first.second ;
                        klass = merge_whatsit( klass, i->second.first.first ) ;
                        is_perfect_match = is_perfect_match && i->second.second ;
                    }

                case 'a':
                    if( verbose >= 2 ) fprintf( stderr, "%s is %s (%d votes)\n", (*s)->id, label[klass], votes ) ;
                    if( verbose >= 3 ) putc( '\n', stderr ) ;
                    summary[klass]++ ;
                    if( is_perfect_match ) ++matched ; else ++mismatched ;
                    break ;

                default:
                    fputs( "don't know how to handle fragment type ", stderr ) ;
                    putc( (*s)->segment, stderr ) ;
                    putc( '\n', stderr ) ;
            }
        }

        {
            double z = 1.96 ; // this is Z_{0.975}, giving a 95% confidence interval (I hope...)
            double k = summary[dirt], n = k + summary[clean] ;
            double p_ = k / n ;
            double c = p_ + 0.5 * z * z / n ;
            double w = z * sqrt( p_ * (1-p_) / n + 0.25 * z * z / (n*n) ) ;
            double d = 1 + z * z / n ;
            double lb = 100.0 * (c-w) / d ;         	// lower bound of CI
            double ml = 100.0 * p_ ;         			// ML estimate
            double ub = 100.0 * (c+w) / d ;      		// upper bound of CI

            for( whatsit klass = unknown ; klass != maxwhatsits ; klass = (whatsit)( (int)klass +1 ) )
            {
                if( mktable ) {
                    printf( "%d\t", summary[klass] ) ;
                } else {
                    printf( "  %s fragments: %d", label[klass], summary[klass] ) ;
                    if( klass == dirt )
                    {
                        printf( " (%.1f .. %.1f .. %.1f%%)", lb, ml, ub ) ;
                    }
                    putchar( '\n' ) ;
                }
            }
            if( mktable ) printf( "%.1f\t%.1f\t%.1f\t", lb, ml, ub ) ;
        }
        { 
            double z = 1.96 ; // this is Z_{0.975}, giving a 95% confidence interval (I hope...)
            double k = mismatched, n = k + matched ;
            double p_ = k / n ;
            double c = p_ + 0.5 * z * z / n ;
            double w = z * sqrt( p_ * (1-p_) / n + 0.25 * z * z / (n*n) ) ;
            double d = 1 + z * z / n ;
            double lb = 100.0 * (c-w) / d ;
            double ml = 100.0 * p_ ;
            double ub = 100.0 * (c+w) / d ;

            if( mktable ) printf( "%d\t%d\t%.1f\t%.1f\t%.1f\n", matched, mismatched, lb, ml, ub ) ;
            else printf( "\n  perfectly aligned fragments: %d"
                         "\n  fragments with differences:  %d (%.1f .. %.1f .. %.1f%%)\n\n",
                         matched, mismatched, lb, ml, ub ) ;
        }
        free_map_alignment( maln ) ;
        free( aln_con ) ;
        free( aln_ass ) ;
    } while(0) ;
}


