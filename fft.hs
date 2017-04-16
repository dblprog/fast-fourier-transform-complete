-- My collaborators are James Waugh, Victor Jiao, Mikhail Iouchkov, Zach Krebs, and Jacob Burroughs.

import Data.Complex

data Poly a = Poly [a]
    deriving (Show)

polyCoeffs :: Poly a -> [a]
polyCoeffs (Poly a) = a

--test for polyCoeffs and for others are at end of file.

polyEval :: RealFloat a => Poly (Complex a) -> (Complex a) -> (Complex a)
polyEval p x = foldl1 ((+) . (*x)) (polyCoeffs p)

polyAddHelper (p) (q)
  | q == [] = p 
  | p == [] = q
  | p == [] && q == [] = []
  | length p < length q = polyAddHelper (0:p) q
  | length q < length p = polyAddHelper p (0:q)
  | otherwise = [(last p) + (last q)] ++ (polyAddHelper ((init p)) ((init q)))

--polyAdd :: RealFloat a => Poly (Complex a) -> (Complex a) -> (Complex a)
polyAdd (p) (q) = Poly (reverse (polyAddHelper (polyCoeffs p) (polyCoeffs q)))

polyScalMult p y = Poly (map (y*) (polyCoeffs p))

polyNegate q = polyScalMult (q) (-1)

polyDiff p q = polyAdd p (polyNegate q)

polyAbs (Poly p) = Poly (map (abs) p)

certError = 0.01
polyApproxHelper p q = abs(p - q) < certError
polyApproxEqual p q = foldr (&&) True $ zipWith (polyApproxHelper) (polyCoeffs p) (polyCoeffs q)

polyLength (Poly p) = length p


equalize xs ys l
  | length xs < l = equalize (0:xs) ys l
  | length ys < l = equalize xs (0:ys) l
  | otherwise = [xs,ys]


setUp2 (Poly a) (Poly b) = equalize (polyCoeffs (Poly a)) (polyCoeffs (Poly b)) (lsum)
  where lsum = polyLength (Poly a) + polyLength (Poly b)

cutDown cs j = take (j+1) (reverse cs)

jthCoeff cs j = foldr (+) 0 $ zipWith (*) (cutDown (cs !! 0) j) (reverse (cutDown (cs !! 1) j))

baseList n = take (n+1) (iterate (+ 1) 0)

polyMult a b = Poly $ reverse $ map (jthCoeff (setUp2 a b)) (baseList (polyLength a + polyLength b - 2))

data PolyPoints x y n = PolyPoints x y n deriving (Show)  
-- Justification: The computational advantage to this order of x and y is that 
-- we can more easily visualize mapping p of x to y.

-- toPolyPoints ::
toPolyPoints p xs n
  | n < polyLength p = undefined
  | otherwise = PolyPoints [xs] (map (polyEval p) [xs]) 

polyPointsAdd (PolyPoints p q n) (PolyPoints r s m)
  | n /= m = undefined
  | otherwise = zipWith (+) q s 

polyPointsScalMult (PolyPoints x y n) c = map (c*) [y]

polyPointsMult (PolyPoints [x] [y] n) (PolyPoints [z] [v] m)
  | [x] /= [z] = undefined
  | otherwise = PolyPoints [x] (zipWith (*) [y] [v]) m

polyPointsApproxEqual x y = foldr (&&) True $ zipWith (polyApproxEqual) (polyCoeffs (Poly x)) (polyCoeffs (Poly y))

evens [] = []
evens [x] = [x]
evens (x:y:zs) = x:(evens zs)

odds [] = []
odds [x] = []
odds (x:y:zs) = y:(odds zs)

data PolyPointsFFT y = PolyPointsFFT [y] 
  deriving (Show)

sumPPFFT l j = sum (fftH l j)

fftH [] j  = []
fftH [a] j = [a]
fftH l j   = (map (ω j) (fftH (odds l) j) ++ (fftH (evens l) j))
  where  
    ω j x = cis (2*pi*(fromIntegral j)/(fromIntegral (length l))) * x

isPow2 l n 
    | length l == n = True
    | length l > n = isPow2 l (n*2)
    | otherwise = False 

fft (Poly l) 
    | isPow2 (polyCoeffs (Poly l)) 2 = PolyPointsFFT $ map (sumPPFFT (pc)) (lst)
    | otherwise = undefined
    where 
          lst = take (polyLength (Poly l)) [0..]
          pc = polyCoeffs (Poly l)

    

invfft (PolyPointsFFT l) = Poly (map (/fromIntegral (length l)) (map (conjugate) (f)))
  where 
    f = getPolyPointsFFTCoeffs (fft (Poly (map (conjugate) l)))
    
getPolyPointsFFTCoeffs :: PolyPointsFFT (Complex a) -> [Complex a]
getPolyPointsFFTCoeffs (PolyPointsFFT a) = a


polyPointsFFTApproxEqual x y = foldr (&&) True $ zipWith (polyApproxEqual) (polyCoeffs (Poly x)) (polyCoeffs (Poly y))


toPolyPointsFFT (Poly p) = fft(polyPad (p))

polyPad p 
  | isPow2 p 2 = (Poly p)
  | otherwise = polyPad (0:p)

polyMultFFT (Poly p) (Poly q) = invfft (PolyPointsFFT (zipWith (*) (g p) (g q)))
  where 
    g w = getPolyPointsFFTCoeffs (g2n (w))
    g2n p = toPolyPointsFFT (Poly (take (maximum [polyLength (Poly p), polyLength (Poly q)]) (repeat 0))) -- getPolyPointsFFTCoeffs (fft (Poly (polyPad p)))
    --a = getPolyPointsFFTCoeffs (toPolyPointsFFT (Poly q)) -- getPolyPointsFFTCoeffs (fft (Poly (polyPad q)))

test1 = polyCoeffs (Poly ([1, 0, 2, 4, 11])) == [1, 0, 2, 4, 11]
test2 = polyEval (Poly [7, 4, 3, 2, 1]) 2 == 161
test3 = polyApproxEqual (polyAdd (Poly [0,1,2,3,4]) (Poly [9, 6, 3])) (Poly[0,1,11,9,7])
test4 = polyApproxEqual(polyScalMult (Poly[3, 4, 5, 6, 13, 2]) 2 ) (Poly[6,8,10,12,26,4])
test5 = polyApproxEqual (polyNegate (Poly [1,2,3,-4])) (Poly [-1,-2,-3,4])
test6 = polyApproxEqual (polyNegate (Poly [1,-2, 3, -4, -5, -7,11 ])) (Poly[-1,2,-3,4,5,7,-11])
test7 = polyApproxEqual (polyDiff (Poly [1,2,3,4,5,6]) (Poly[0,1,2,3,4,5])) (Poly [1,1,1,1,1,1])
test8 = polyApproxEqual (polyMult (Poly [4,3,1]) (Poly [1,2]))  (Poly[4,11,7,2])
{-
test9 = polyPointsApproxEqual (toPolyPoints [1,1,1,1] [1,2,3,4] 4) (PolyPoints [1, 2, 3, 4] [4 15 40 85] 4)
test10 = polyApproxEqual (Poly (evens [1,2,3,4,5,6,7,8])) (Poly (odds [1,2,3,4,5,6,7,8,9]))
test11 = polyApproxEqual (Poly odds [6,5,4,3,2,1]) (Poly [6,4,2])
test12 = polyApproxEqual (Poly evens [6,5,4,3,2,1]) (Poly [5,3,1]) 
test13 = polyPointsFFTApproxEqual (toPolyPointsFFT (Poly [1,2,3,4])) fft(polyPad [1,2,3,4])
test14 = polyApproxEqual (Poly (invfft (fft (Poly [1,2,3,4])))) (Poly [1,2,3,4])
-}
tests = test1 && test2 && test3 && test4 && test5 && test6 && test7 && test8 -- && test9 && test10 && test11 && test12 && test13 && test14

