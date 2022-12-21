{-# LANGUAGE DeriveGeneric #-}
module Main (main) where

import ParBnC (
    ObjectiveType (Maximization, Minimization),
    BranchProblem,
    toTableau,
    addSlackMatrix,
    constructBranchAndBound,
    constructParBranchAndBound,
    constructBranchAndCut,
    constructParBranchAndCut,
    searchBBTreeMax)
import Data.Aeson
    ( (.:),
      object,
      FromJSON(parseJSON),
      Value(Object),
      KeyValue((.=)),
      ToJSON(toJSON),
      eitherDecode
    )
import Data.Text (Text)
import qualified Data.ByteString.Lazy as B
import GHC.Generics
import System.Exit (die)
import System.Environment (getArgs, getProgName)
import Numeric.LinearAlgebra as LA
    ( fromLists, fromList, Matrix, R, Vector, size)

data IP = IP {
    name :: !Text,
    costC :: [Double],
    matA :: [[Double]],
    boundB :: [Double]
} deriving (Show, Generic)

instance FromJSON IP 
instance ToJSON IP

jsonFile :: FilePath
jsonFile = "test/medium_ip.json"

getJSON :: IO B.ByteString
getJSON = B.readFile jsonFile

obj :: ObjectiveType
obj = Maximization

main :: IO ()
main = do
    args <- getArgs
    case args of
        -- Sequential
        [fileName, maxDepth, "seq", bc, to] -> do
            let jsonFile = fileName
            let getJSON = B.readFile jsonFile
            let md = read maxDepth :: Int
            d <- (eitherDecode <$> getJSON) :: IO (Either String IP)
            case d of
                Left err -> putStrLn $ "Failed to load test case" ++ err
                Right ps -> do
                    let testC = LA.fromList $ costC ps
                    let testMat = LA.fromLists $ matA ps
                    let testB = LA.fromList $ boundB ps
                    let testTab = addSlackMatrix $ toTableau testC testMat testB
                    let testIntMask = LA.fromList $ replicate (LA.size testC) True
                    case bc of
                        -- Branch and Bound
                        "b" -> do
                            case to of
                                -- full tree output
                                "tree" -> print $ constructBranchAndBound obj testTab testIntMask testC md
                                -- just optimal solution
                                "solution" -> print $ searchBBTreeMax $ constructBranchAndBound obj testTab testIntMask testC md
                                _ -> error "Wrong input format for output type"
                        -- Branch and Cut
                        "c" -> do 
                            case to of
                                -- full tree output
                                "tree" -> print $ constructBranchAndCut obj testTab testIntMask testC md
                                -- just optimal solution
                                "solution" -> print $ searchBBTreeMax $ constructBranchAndCut obj testTab testIntMask testC md
                        _ -> error "Wrong input format for Bound/Cut" 
        [fileName, maxDepth, "par", bc, to, parDepth] -> do
            let jsonFile = fileName
            let getJSON = B.readFile jsonFile
            let md = read maxDepth :: Int
            let pd = read parDepth :: Int
            d <- (eitherDecode <$> getJSON) :: IO (Either String IP)
            case d of
                Left err -> putStrLn $ "Failed to load test case" ++ err
                Right ps -> do
                    let testC = LA.fromList $ costC ps
                    let testMat = LA.fromLists $ matA ps
                    let testB = LA.fromList $ boundB ps
                    let testTab = addSlackMatrix $ toTableau testC testMat testB
                    let testIntMask = LA.fromList $ replicate (LA.size testC) True
                    case bc of
                        -- Branch and Bound
                        "b" -> do
                            case to of
                                -- full tree output
                                "tree" -> print $ constructParBranchAndBound obj testTab testIntMask testC md pd
                                -- just optimal solution
                                "solution" -> print $ searchBBTreeMax $ constructParBranchAndBound obj testTab testIntMask testC md pd
                                _ -> error "Wrong input format for output type"
                        -- Branch and Cut
                        "c" -> do 
                            case to of
                                -- full tree output
                                "tree" -> print $ constructParBranchAndCut obj testTab testIntMask testC md pd
                                -- just optimal solution
                                "solution" -> print $ searchBBTreeMax $ constructParBranchAndCut obj testTab testIntMask testC md pd
                        _ -> error "Wrong input format for Bound/Cut" 
        _ -> do
            pn <- getProgName
            die $ "Usage: stack exec " ++ pn 
                ++ "<json_filename> <maxDepth> <seq> <bc> <to> <parDepth>"
