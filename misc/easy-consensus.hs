{-# LANGUAGE BangPatterns, NamedFieldPuns, RecordWildCards #-}
import Control.Monad
import Data.Char ( toUpper )
import Data.List ( transpose, sortBy, sort )
import System.Console.GetOpt
import System.Environment ( getArgs, getProgName )
import System.Exit
import System.IO

data Conf = Opts { percent :: Double, to_ambicode :: String -> Char, output :: String -> IO () }

defaultOpts :: IO Conf
defaultOpts = return $ Opts 1 toAmbicode putStr

options :: [OptDescr (Conf -> IO Conf)]
options = [
    Option "p" ["percent"] (ReqArg (\a c -> readIO a >>= \p -> return $ c { percent = p / 100 }) "P")
                           "Set percentage needed for consensus to P",
    Option "n" ["only-n"]  (NoArg (\c -> return $ c { to_ambicode = toNucleotide }))
                           "Allow only nucleotides and 'N' in consensus",
    Option "i" ["iupac"]   (NoArg (\c -> return $ c { to_ambicode = toUpper . toAmbicode }))
                           "Allow all IUPAC ambiguity codes in consensus",
    Option "g" ["gaps"]    (NoArg (\c -> return $ c { to_ambicode = toAmbicode }))
                           "Allow IUPAC codes and small letters for optional gaps",
    Option "o" ["output"]  (ReqArg (\fn c -> return $ c { output = if fn == "-" then putStr else writeFile fn }) "FILE")
                           "Write output to FILE instead of stdout",
    Option "h?" ["help","usage"] (NoArg usage)
                           "Print this helpful message" ]
  where
    usage _ = do pn <- getProgName
                 let blah = "Usage: " ++ pn ++ "[options...] [fasta-file]\n\
                            \Reads a multi-FastA file, computes a consensus where \
                            \a given fraction of the sequences agree, \
                            \and writes it out in FastA format."
                 hPutStrLn stderr $ usageInfo blah options
                 exitSuccess

main :: IO ()
main = do
    (opts', files, errors) <- getOpt Permute options `fmap` getArgs
    unless (null errors) $ mapM_ (hPutStrLn stderr) errors >> exitFailure
    Opts{..} <- foldl (>>=) defaultOpts opts'

    output . writeFasta . filter (/= '-') .
        map (call_cons percent to_ambicode . map toUpper) . 
        transpose . map snd . concat =<<
        mapM (\fn -> readMFasta `fmap` if fn == "-" then getContents else readFile fn) files


readMFasta :: String -> [(String, String)]
readMFasta = go . lines . filter (/= '\r')
  where
    go ls = case dropWhile (not . isHeader) ls of
                [] -> []
                (h:ls1) -> let (b,ls2) = span (not.isHeader) ls1
                           in (h, concat b) : go ls2

writeFasta :: String -> String
writeFasta s = unlines $ ">consensus" : split 60 s
  where
    split n [] = []
    split n s' = case splitAt n s' of (l,r) -> l : split n r

isHeader :: String -> Bool
isHeader ('>':_) = True
isHeader _ = False

-- Consensus call into ambiguity code:
-- - count A,C,G,T,-; sort by prevalence
-- - take the more common stuff until we reach a fraction of p
-- - turn into ambiguity code, small case if the - is included
call_cons :: Double -> (String -> Char) -> String -> Char
call_cons p to_ambicode s 
    | null cum = 'N'
    | otherwise = case span ((< ceiling (p * total)) . fst) cum of
                    (l,r:_) -> collapse (r:l)
                    (l,[ ]) -> collapse l
  where
    total = fromIntegral $ fst $ last cum
  
    cum = scanl1 (\(a,_) (b,n) -> (a+b,n)) $
          sortBy (\(a,_) (b,_) -> compare b a) $
          zip [a,c,g,t,z] "ACGT-"

    (a,c,g,t,z,n) = count 0 0 0 0 0 0 s

    collapse = to_ambicode . sort . map snd

toAmbicode :: String -> Char
toAmbicode ""     = '-'
toAmbicode "A"    = 'A'
toAmbicode "AC"   = 'M'
toAmbicode "ACG"  = 'V'
toAmbicode "ACGT" = 'N'
toAmbicode "ACT"  = 'G'
toAmbicode "AG"   = 'R'
toAmbicode "AGT"  = 'D'
toAmbicode "AT"   = 'W'
toAmbicode "C"    = 'C'
toAmbicode "CG"   = 'S'
toAmbicode "CGT"  = 'B'
toAmbicode "CT"   = 'Y'
toAmbicode "G"    = 'G'
toAmbicode "GT"   = 'K'
toAmbicode "T"    = 'T'

toAmbicode "-"     = '-'
toAmbicode "-A"    = 'a'
toAmbicode "-AC"   = 'm'
toAmbicode "-ACG"  = 'v'
toAmbicode "-ACGT" = 'n'
toAmbicode "-ACT"  = 'h'
toAmbicode "-AG"   = 'r'
toAmbicode "-AGT"  = 'd'
toAmbicode "-AT"   = 'w'
toAmbicode "-C"    = 'c'
toAmbicode "-CG"   = 's'
toAmbicode "-CGT"  = 'b'
toAmbicode "-CT"   = 'y'
toAmbicode "-G"    = 'g'
toAmbicode "-GT"   = 'k'
toAmbicode "-T"    = 't'

toAmbicode x = error $ "huh?  " ++ show x

toNucleotide :: String -> Char
toNucleotide [c] | c `elem` "ACGT-" = c
toNucleotide [ ] = '-'
toNucleotide  _  = 'N'


count :: Int -> Int -> Int -> Int -> Int -> Int
      -> String -> (Int,Int,Int,Int,Int,Int)
count !a !c !g !t !z !n s = case s of
    [      ] -> (a,c,g,t,z,n)
    ('A':s') -> count (a+1) c g t z (n+1)  s'
    ('C':s') -> count a (c+1) g t z (n+1)  s'
    ('G':s') -> count a c (g+1) t z (n+1)  s'
    ('T':s') -> count a c g (t+1) z (n+1)  s'
    ('-':s') -> count a c g t (z+1) (n+1)  s'
    ( _ :s') -> count a c g t z     (n+1)  s'

