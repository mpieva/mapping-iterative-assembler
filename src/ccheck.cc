#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <utility>
#include <getopt.h>

extern "C" {
#include "map_align.h"
#include "mia.h"
}

#include "myers_align.h"

/*
 * Rough cut at contamination checker.  Outline:
 *
 * - read the human reference; this will contain ambiguity codes
 * - read maln file, including assembly and assembled reads
 * - align reference-consensus and assembly globally
 *   This uses my O(nd) aligner, for it grasps ambiguity codes and runs
 *   fast enough for long, but similar sequences.
 * - find "diagnostic positions", positions where ass and ref differ
 * - for every "end" fragment: store it (or its summary?) and
 *   later join with its other half
 * - for every "full" fragment: if it crosses at least one diagnostic
 *   position, cut out that range from ref and align to it globally
 * - for every position where the bases agree, classify it, then
 *   classify the fragment (conflicting, uninformative, contaminant,
 *   endogenous)
 * - produce a summary
 *
 * Additional features:
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
 *   everywhere sensibly, then use it
 * - somehow account for sequencing error, especially on Solexa
 * - remove linear scans where binary search methods would work
 */

namespace std {
template <typename U, typename V>
std::ostream& operator << ( std::ostream& s, const std::pair<U,V> p )
{
	return s << '<' << p.first << ',' << p.second << '>' ;
}
}

void print_aln( const char* aln )
{
	while( *aln ) {
		int p = 0 ;
		for( const char *a = aln ; *a && p != 72 ; ++a, ++a, ++p )
			std::cout << *a ;
		std::cout << '\n' ;
		p = 0 ;
		for( ; *aln && p != 72 ; ++aln, ++aln, ++p )
			std::cout << aln[1] ;
		std::cout << '\n' ;
		std::cout << '\n' ;
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
bool is_diagnostic( const char* aln )
{
	return aln[0] != aln[1]
				&& aln[0] != 'N' && aln[1] != 'N' 
				&& aln[0] != '-' && aln[1] != '-' ;
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
bool is_transversion_( std::pair< int, std::pair< char, char > > p ) 
{ return is_transversion( p.second.first, p.second.second ) ; }


dp_list mk_dp_list( const char* aln, bool transversions, int span_from, int span_to )
{
	dp_list l ;
	for( int p = span_from ; p != span_to && *aln ; ++aln, ++aln )
	{
		if( is_diagnostic( aln ) && ( !transversions || is_transversion( aln[0], aln[1] )))
			l[p] = std::make_pair( aln[0], aln[1] ) ;
		if( aln[1] != '-' ) ++p ;
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
std::string lift_over( const char* aln, int s, int e ) 
{
	std::string r ;
	for( int p = 0 ; p < e && *aln ; ++aln, ++aln )
	{
		if( aln[0] != '-' && p >= s ) r.push_back( aln[0] ) ;
		if( aln[1] != '-' ) ++p ;
	}
	return r ;
}
 
std::string compare_cstr_str( const char* pc, const std::string& s ) 
{
	std::string r ;
	for( std::string::const_iterator b = s.begin(), e = s.end() ;
			b != e && *pc ; ++b, ++pc )
		r.push_back( *b == *pc ? '*' : ' ' ) ;
	return r ;
}

bool consistent( bool adna, char x, char y )
{
	char x_ = x == 'G' ? 'R' : x == 'C' ? 'Y' : x ;
	return x == '-' || y == '-' || (char_to_bitmap( adna ? x_ : x ) & char_to_bitmap(y)) != 0 ;
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

struct option longopts[] = {
	{ "reference", required_argument, 0, 'r' },
	{ "ancient", no_argument, 0, 'a' },
	{ "verbose", no_argument, 0, 'v' },
	{ "help",    no_argument, 0, 'h' },
	{ "transversions", no_argument, 0, 't' },
	{ "span", required_argument, 0, 's' },
	{ "maxd", required_argument, 0, 'd' },
	{ 0,0,0,0 }
} ;

void usage( const char* pname )
{
	std::cout << "Usage: " << pname << " [-r <ref.fa>] [-a] [-v] <aln.maln> \n\n"
		"Reads a maln file and tries to quantify contained contamination.\n"
		"Options:\n"
		"  -r, --reference FILE     FASTA file with the likely contaminant\n"
		"  -a, --ancient            treat DNA as ancient (i.e. likely deaminated)\n"
		"  -t, --transversions      only transversions are diagnostic\n"
		"  -s, --span M-N           only look at range from M to N\n"
		"  -d, --maxd D             allow up to D differences between the references\n"
		"  -v, --verbose            increases verbosity level (can be repeated)\n"
		"  -h, --help               print this help message\n" << std::endl ;
}

int main( int argc, char * const argv[] )
{
	int summary[ maxwhatsits ] = {0} ;
	struct refseq hum_ref ;
	bool adna = false ;
	bool transversions = false ;
	int verbose = 0 ;
	int maxd = 1000 ;
	const char* ref_file = "mt311.fna" ;
	int span_from = 0, span_to = INT_MAX ;

	if( argc == 0 ) { usage( argv[0] ) ; return 0 ; }

	int opt ;
	do {
		opt = getopt_long( argc, argv, "r:avhts:d:", longopts, 0 ) ;
		switch( opt ) 
		{
			case 'r': 
				ref_file = optarg ;
				break ;
			case 'a':
				adna = true ;
				break ;
			case 'v':
				++verbose ;
				break ;
			case ':':
				std::clog << "missing option argument" << std::endl ;
				break ;
			case '?':
				std::clog << "unknown option" << std::endl ;
				break ;
			case 'h':
				usage( argv[0] ) ;
				return 0 ;
			case 't':
				transversions = true ;
				break ;
			case 's':
				sscanf( optarg, "%u-%u", &span_from, &span_to ) ;
				span_from-- ;
				break ;
			case 'd':
				maxd = atoi( optarg ) ;
				break ;
		}
	} while( opt != -1 ) ;

	if( optind == argc ) { usage( argv[0] ) ; return 0 ; }

	read_fasta_ref( &hum_ref, ref_file ) ;
	MapAlignmentP maln = read_ma( argv[optind] ) ;
	PSSMP submat = maln->fpsm ;

	char aln[ ( strlen(hum_ref.seq) + strlen(maln->ref->seq) + 1 ) *2 ] ;
	int d = myers_align( hum_ref.seq, maln->ref->seq, maxd, equal, aln ) ;

	if( d == maxd ) { std::cout << "Couldn't align references (try to increase maxd)." << std::endl ; return 1 ; }
	if( verbose >= 1 ) std::cout << d << " total differences between reference and assembly." << std::endl ;
	if( verbose >= 6 ) print_aln( aln ) ;
	
	dp_list l = mk_dp_list( aln, transversions, span_from, span_to ) ;

	if( verbose >=1 ) 
		std::cout << l.size() << " diagnostic positions, "
		<< std::count_if( l.begin(), l.end(), is_transversion_ ) 
		<< " of which are transversions.\n" ;
	if( verbose >= 3 ) 
	{
		std::copy( l.begin(), l.end(), std::ostream_iterator
				<dp_list::value_type>( std::cout, ", " ) ) ;
		std::cout << std::endl ;
	}

	typedef std::map< std::string, std::pair< whatsit, int > > Bfrags ;
	Bfrags bfrags ;

	for( const AlnSeqP *s = maln->AlnSeqArray ;
			s != maln->AlnSeqArray + maln->num_aln_seqs ; ++s )
	{
		whatsit klass = unknown ;
		int votes = 0 ;

		std::string seq_id = (*s)->id ;
		seq_id.push_back('/') ;
		seq_id.push_back( (*s)->segment ) ;

		std::pair< dp_list::const_iterator, dp_list::const_iterator > p =
			overlapped_diagnostic_positions( l, *s ) ;
		if( p.first == p.second ) {
			if( verbose >= 3 ) std::cout << seq_id << ": no diagnostic positions" << std::endl ;
		}
		else
		{
			if( verbose >= 3 )
			{
				std::cout << seq_id << ": " << std::distance( p.first, p.second )
				          << " diagnostic positions" ;
				if( verbose >= 4 ) 
				{
					std::cout << ": " ;
					std::copy( p.first, p.second, std::ostream_iterator
							< dp_list::value_type >( std::cout, " " ) ) ;
				}
				std::cout << std::endl ;
				std::cout << "range:  " << (*s)->start << ".." << (*s)->end << std::endl ;
			}

			std::string the_read ;
			for( char *nt = (*s)->seq, **ins = (*s)->ins ; *nt ; ++nt, ++ins )
			{
				if( *nt != '-' ) the_read.push_back( *nt ) ;
				if( *ins ) the_read.append( *ins ) ;
			}
			std::string the_ass( maln->ref->seq + (*s)->start, (*s)->end - (*s)->start + 1 ) ;
			std::string lifted = lift_over( aln, (*s)->start, (*s)->end + 1 ) ;

			if( verbose >= 5 )
				std::cout << "raw read: " << the_read << '\n'
					<< "lifted:   " << lifted << '\n'
					<< "assembly: " << the_ass << '\n'
					<< '\n'
					<< "aln.read: " << (*s)->seq << '\n'
					<< "aln.assm: " << the_ass << '\n' 
					<< "matches:  " << compare_cstr_str( (*s)->seq, the_ass ) << '\n' ;

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
				std::cout << '\n'
					<< "aln.read: " << pwaln.frag_seq << '\n'
					<< "aln.ref:  " << pwaln.ref_seq << '\n' 
					<< "matches:  " << compare_cstr_str( pwaln.frag_seq, std::string(pwaln.ref_seq) ) << "\n\n" ;

			free_alignment( frag_aln ) ;

			char *paln = aln ;
			int ass_pos = 0 ;
			while( ass_pos != (*s)->start && *paln ) 
			{
				if( paln[1] != '-' ) ass_pos++ ;
				paln += 2 ;
			}

			std::string in_ref = lifted.substr( 0, pwaln.start ) ;
			in_ref.append( pwaln.ref_seq ) ;

			char *in_frag_v_ref = pwaln.frag_seq ;
			char *in_ass = maln->ref->seq + (*s)->start ;
			char *in_frag_v_ass = (*s)->seq ;

			if( paln[0] != in_ref[0] || paln[0] == '-' )
			{
				std::cout << "huh? (R+" << pwaln.start << ") " ;
				for( const char* pc = paln ; *pc && pc != paln+20 ; pc+=2 ) std::cout << *pc ;
				std::cout << ' ' << in_ref.substr(0,10) << std::endl ;
			}
			if( paln[1] != in_ass[0] && paln[1] != '-' )
			{
				std::cout << "huh? (A+" << pwaln.start << ") " ;
				for( const char* pc = paln+1 ; *pc && pc != paln+21 ; pc+=2 ) std::cout << *pc ;
				std::cout << ' ' << std::string( in_ass, in_ass+10 ) << std::endl ;
			}

			while( ass_pos != (*s)->end +1 && *paln && !in_ref.empty() && *in_ass && *in_frag_v_ass && *in_frag_v_ref )
			{
				if( is_diagnostic( paln ) ) {
					if( verbose >= 4 )
						std::cout << "diagnostic pos.: " << ass_pos << ' ' 
							<< in_ref[0] << '/' << *in_frag_v_ref << ' '
							<< *in_ass << '/' << *in_frag_v_ass << ' ' ;
					if( *in_frag_v_ref != *in_frag_v_ass ) 
					{
						if( verbose >= 4 ) std::cout << "in disagreement." << std::endl ;
					}
					else
					{
						bool maybe_clean = consistent( adna, *in_ass, *in_frag_v_ass ) ;
						bool maybe_dirt =  consistent( adna, in_ref[0], *in_frag_v_ref ) ;

						if( verbose >= 4 )
							std::cout << ( maybe_dirt  ? "" : "in" ) << "consistent/"
								<< ( maybe_clean ? "" : "in" ) << "consistent" 
								<< std::endl ;

						if( maybe_clean && !maybe_dirt && klass == unknown ) klass = clean ;
						if( maybe_clean && !maybe_dirt && klass == dirt    ) klass = conflict ;
						if( !maybe_clean && maybe_dirt && klass == unknown ) klass = dirt ;
						if( !maybe_clean && maybe_dirt && klass == clean   ) klass = conflict ;
						if( !maybe_clean && !maybe_dirt )                    klass = nonsense ;
						if( maybe_dirt != maybe_clean ) votes++ ;
					}
				}

				if( paln[0] != '-' ) {
					do {
						in_ref=in_ref.substr(1) ;
						in_frag_v_ref++ ;
					} while( in_ref[0] == '-' ) ;
				}
				if( paln[1] != '-' ) {
					ass_pos++ ;
					do {
						in_ass++ ;
						in_frag_v_ass++ ;
					} while( *in_ass == '-' ) ;
				}
				paln += 2 ;
			}
		}

		Bfrags::const_iterator i = bfrags.find( (*s)->id ) ;

		switch( (*s)->segment )
		{
			case 'b':
				bfrags[ (*s)->id ] = std::make_pair( klass, votes ) ;
				break ;

			case 'f':
				if( i == bfrags.end() ) 
				{
					std::clog << seq_id << " is missing its back." << std::endl ;
				}
				else
				{
					votes += i->second.second ;
					klass = merge_whatsit( klass, i->second.first ) ;
				}
				
			case 'a':
				if( verbose >= 2 ) 
					std::cout << seq_id << " is " << label[klass]
						<< " (" << votes << " votes)" << std::endl ;
				if( verbose >= 3 ) std::cout << std::endl ;
				summary[klass]++ ;
				break ;

			default:
				std::clog << "don't know how to handle fragment type " << (*s)->segment << std::endl ;
		}
	}

	std::cout << "\nSummary: \n" ;
	for( whatsit klass = unknown ; klass != maxwhatsits ; klass = (whatsit)( (int)klass +1 ) )
	{
		std::cout << label[klass] << " fragments: " << summary[klass] ;
		if( klass == dirt )
		{
			double z = 1.96 ; // this is Z_{0.975}, giving a 95% confidence interval (I hope...)
			double k = summary[dirt], n = k + summary[clean] ;
			double p_ = k / n ;
			double c = p_ + 0.5 * z * z / n ;
			double w = z * std::sqrt( p_ * (1-p_) / n + 0.25 * z * z / (n*n) ) ;
			double d = 1 + z * z / n ;

			std::cout << " (" << std::setprecision(1)
				<< 100.0 * (c-w) / d << " .. "		// lower bound of CI
				<< 100.0 * p_ << " .. "				// ML estimate
				<< 100.0 * (c+w) / d << "%)" ;		// upper bound of CI
		}
		std::cout << '\n' ;
	}
	std::cout << std::endl ;
}


