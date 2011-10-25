-------------------------------------------------------------------------------
--                                                                           --
--                              PI Calculation                               --
--                                                                           --
--                            PI_Calculation.ads                             --
--                                                                           --
--                                  SPEC                                     --
--                                                                           --
--                   Copyright (C) 1996 Ulrik Hørlyk Hjort                   --
--                                                                           --
--  PI Calculation is free software;  you can  redistribute it               --
--  and/or modify it under terms of the  GNU General Public License          --
--  as published  by the Free Software  Foundation;  either version 2,       --
--  or (at your option) any later version.                                   --
--  PI Calculation is distributed in the hope that it will be                --
--  useful, but WITHOUT ANY WARRANTY;  without even the  implied warranty    --
--  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  --
--  See the GNU General Public License for  more details.                    --
--  You should have  received  a copy of the GNU General                     --
--  Public License  distributed with Yolk.  If not, write  to  the  Free     --
--  Software Foundation,  51  Franklin  Street,  Fifth  Floor, Boston,       --
--  MA 02110 - 1301, USA.                                                    --
--                                                                           --
-------------------------------------------------------------------------------
package PI_Calculation is

   type Point_T is
      record
         X : Long_Long_Float;
         Y : Long_Long_Float;
      end record;

   --------------------------------------------------------------------------------------
   -- PI Approximated by Monte Carlo method given by:
   --
   -- A circle C is given by center (0,0) and radius R is inscribed within a square S with
   -- side length 2*R.
   --
   -- Area for C is: AC = PI * (R ** 2)
   --
   -- Area for S is: SC = (2 * R) ** 2
   --
   -- The ratio Rho of AC to SC is: Rho = AC/AS = PI * (R**2) / 4 (R**2) = PI/4
   --
   -- Some points in the unit square are random selected where M is the number of
   -- points inside the unit circle satisfying X**2+Y**2 <= 1, and N is the number of
   -- random points picked.
   -- Then we get Rho = M/N and can approximate PI by Rho * 4.
   --------------------------------------------------------------------------------------
   procedure Monte_Carlo_Method(Epsilon : in Long_Long_Float);

   --------------------------------------------------------------------------------------
   -- PI  approximation bt the Gregory Leibniz Series:
   --
   --  4 * (1/1 - 1/3 + 1/5 - 1/7 + 1/9 - 1/11 ...)
   --
   --------------------------------------------------------------------------------------
   procedure Gregory_Leibniz_Series(Epsilon : in Long_Long_Float);

   --------------------------------------------------------------------------------------
   -- PI Approximation by Wallis product formula:
   --
   -- PI/2 = (2/1) * (2/3) * (4/3) * (4/5) ... = 2 * PROD((N**2) / ((N**2) -1)),
   --                                           for N mod 2 = 0
   --
   -- Source: The American Mathematical Monthly November 1993
   --------------------------------------------------------------------------------------
   procedure Wallis_Product(Epsilon : in Long_Long_Float);

   --------------------------------------------------------------------------------------
   --  Binary PI digit computation by the Bailey Borwein Plouffe formula:
   --
   --  Pi = SUM ( 1/(16 ** K) * ( 4/(8*K + 1) - 2/(8*K + 4) - 1/(8*K + 5) - 1/(8*K + 6))),
   --        Where K in [0, inf]
   --
   --------------------------------------------------------------------------------------
   procedure Bailey_Borwein_Plouffe(Digit : Natural);

   --------------------------------------------------------------------------------------
   --
   --  PI digit computation variant by the Simon Plouffe
   --
   --------------------------------------------------------------------------------------
   function Plouffe(Digit : Natural) return Long_Long_Integer;
end PI_Calculation;
