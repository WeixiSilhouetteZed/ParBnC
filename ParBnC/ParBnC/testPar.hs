import Control.Parallel(par)
import GHC.Integer (Integer)

nfib2 :: Integer -> Integer 
nfib2 n | n < 2 = 1
nfib2 n = par nf (nf + nfib2 (n - 2) + 1)
    where nf = nfib2 (n - 1)

main :: IO()
main = print(nfib2 40)