-------------------------------------------------------------------------------
--                                                                           --
--                              PI Calculation                               --
--                                                                           --
--                            PI_Calculation.adb                             --
--                                                                           --
--                                  BODY                                     --
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
with Ada.Numerics.Discrete_Random; use Ada.Numerics;
with Ada.Text_IO; use Ada.Text_IO;
with Ada.Float_Text_IO; use Ada.Float_Text_IO;
with Ada.Numerics.Elementary_Functions; use Ada.Numerics.Elementary_Functions;
with Number_Theory_Tools; use Number_Theory_Tools;

package body PI_Calculation is

   Pi_C : constant Long_Long_Float := 3.14159265358979323846;

   subtype Random_Interval is Natural range 0..Natural'Last;
   package Random_Natural is new Discrete_Random(Random_Interval);
   use Random_Natural;

   G : Generator;

   --------------------------------------------------------------------------------------
   --
   -- Returns true if Point inside Circle with Radius and center (0,0)
   --
   --
   --------------------------------------------------------------------------------------
   function Point_Inside_Circle(Point : Point_T; Radius : Long_Long_Float) return Boolean is
   begin
      if (Point.X ** 2) + (Point.Y ** 2) < (Radius ** 2) then
         return True;
      else
         return False;
      end if;
   end Point_Inside_Circle;


   --------------------------------------------------------------------------------------
   --
   -- Returns a random point (X,Y) with  0 =< X <= Upper_limit /\ 0 =< Y <= Upper_limit
   --
   --------------------------------------------------------------------------------------
   function Random_Point(Upper_Limit : Long_Long_Float) return Point_T is

      Random_Point    : Point_T;

   begin
      Random_Point.X := Long_Long_Float(Random(G))/Long_Long_Float(Natural'Last) * Upper_Limit;
      Random_Point.Y := Long_Long_Float(Random(G))/Long_Long_Float(Natural'Last) * Upper_Limit;

      return Random_Point;
   end Random_Point;



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
   procedure Monte_Carlo_Method(Epsilon : in Long_Long_Float) is
      Radius : constant Long_Long_Float := 1.0;
      --      Number_Of_Points : constant Positive := Positive'Last;
      Number_Of_Points : constant Positive := 9999999;
      Points_Inside_Circle : Natural := 0;
      Point : Point_T;
      PI_Approx : Long_Long_Float := 0.0;

   begin
      for N in 1 .. Number_Of_Points Loop
         Point := Random_Point(Radius);
         if Point_Inside_Circle(Point, Radius) then
            Points_Inside_Circle := Points_Inside_Circle + 1;
         end if;
         PI_Approx := (4.0 * Long_Long_Float(Points_Inside_Circle)) / Long_Long_Float(N);
         exit when (abs((4.0 * PI_Approx) - PI_C) < Epsilon);

      end loop;

      Put("PI Approx by Monte Carlo method = ");
      Put(Item => Float(PI_Approx), Fore => 1, Aft => 18, EXP => 0);
      Put( " with Error = " & Float'Image(Float(abs((4.0 * PI_Approx) - PI_C))));
      New_Line;
   end Monte_Carlo_Method;


   --------------------------------------------------------------------------------------
   -- PI  approximation bt the Gregory Leibniz Series:
   --
   --  4 * (1/1 - 1/3 + 1/5 - 1/7 + 1/9 - 1/11 ...)
   --
   --------------------------------------------------------------------------------------
   procedure Gregory_Leibniz_Series(Epsilon : in Long_Long_Float) is
      PI_Approx : Long_Long_Float := 0.0;
      N   : Long_Long_Float := 1.0;
   begin
      loop
         PI_Approx := PI_Approx + ((-1.0 *((-1.0) ** Natural(N))) * (1.0/((2.0* N)-1.0)));
         exit when (abs((4.0 * PI_Approx) - PI_C) < Epsilon);
         N := N + 1.0;
      end loop;
      PI_Approx := PI_Approx * 4.0;

      Put("PI Approx by Gregory-Leibniz Series = ");
      Put(Item => Float(PI_Approx), Fore => 1, Aft => 18, EXP => 0);
      Put(" after" & Integer'Image(Integer(N)) & " Iterations");
      Put( " with Epsilon = " & Float'Image(Float(Epsilon)));
      New_Line;
   end Gregory_Leibniz_Series;


   --------------------------------------------------------------------------------------
   -- PI Approximation by Wallis product formula:
   --
   -- PI/2 = (2/1) * (2/3) * (4/3) * (4/5) ... = 2 * PROD((N**2) / ((N**2) -1)),
   --                                           for N mod 2 = 0
   --
   -- Source: The American Mathematical Monthly November 1993
   --------------------------------------------------------------------------------------
   procedure Wallis_Product(Epsilon : in Long_Long_Float) is
      PI_Approx : Long_Long_Float := 1.0;
      N   : Long_Long_Float := 2.0;
   begin
      loop
         PI_Approx := PI_Approx * ((N ** 2) / ((N ** 2) - 1.0));
         exit when (abs((2.0 * PI_Approx) - PI_C) < Epsilon);
         --      Put(Item => Float(PI_Approx *2.0), Fore => 1, Aft => 18, EXP => 0);New_Line;
         N := N + 2.0;
      end loop;
      PI_Approx := PI_Approx * 2.0;

      Put("PI Approx by Wallies Product = ");
      Put(Item => Float(PI_Approx), Fore => 1, Aft => 18, EXP => 0);
      Put(" after" & Integer'Image(Integer(N)) & " Iterations");
      Put( " with Epsilon = " & Float'Image(Float(Epsilon)));
      New_Line;
   end Wallis_Product;


   --------------------------------------------------------------------------------------
   --
   --
   --
   --
   --------------------------------------------------------------------------------------
   function Exponentiation_Scheme(P : Long_Long_Float; Ak : Long_Long_Float) return Long_Long_Float is
      Limit : constant Natural := 24;
      type Float_Array_T is array (Natural range <>) of Long_Long_Float;
      Powers : Float_Array_T(0..Limit) := (others => 0.0);
      K      : Natural := 1;
      P1     : Long_Long_Float;
      P2     : Long_Long_Float;
      R      : Long_Long_Float := 1.0;

   begin
      if Ak = 1.0 then
         return 0.0;
      else
         Powers(0) := 1.0;
         for I in 1 .. Limit-1 loop
            Powers(I) := 2.0 * Powers(I-1);
         end loop;

         loop
            exit when Powers(K) > P;
            K := K + 1;
         end loop;

         P2 := Powers(K-1);
         P1 := P;

         loop
            exit when P2 < 1.0;
            if P1 >= P2 then
               R := F_Mod(16.0 * R, Ak);
               P1 := P1 - P2;
            end if;
            P2 := P2 * 0.5;
            if P2 >= 1.0 then
               R := F_Mod(R*R,Ak);
            end if;
         end loop;
         return F_Mod(R,Ak);
      end if;
   end Exponentiation_Scheme;

   --------------------------------------------------------------------------------------
   --
   --
   --
   --
   --------------------------------------------------------------------------------------
   function Evaluate_Serie(N : Long_Long_Float; Digit : Natural) return Long_Long_Float is
      Ak : Long_Long_Float := 0.0;
      P  : Long_Long_Float := 0.0;
      S  : Long_Long_Float := 0.0;
      T  : Long_Long_Float := 0.0;

      Digit_Tmp : Long_Long_Float := Long_Long_Float(Digit);

      Epsilon : constant Long_Long_Float := 1.0E-17;

   begin
      for K in 0 .. Digit - 1 loop
         Ak := 8.0 * Long_Long_Float(K) + N;
         P  := Long_Long_Float(Digit - K);
         T  := Exponentiation_Scheme(P,Ak);
         S  := S + (T/Ak);
         S  := S - Long_Long_Float'Floor(S);
      end loop;

      P := 16.0;

      loop
         Ak := 8.0 * Digit_Tmp + N;
         P := P / 16.0;
         T := P / Ak;
         S := S + T;
         S  := S - Long_Long_Float'Floor(S);
         exit when T <= Epsilon;
         Digit_Tmp := Digit_Tmp + 1.0;
      end loop;
      return S;
   end Evaluate_Serie;


   --------------------------------------------------------------------------------------
   --
   --
   --
   --
   --------------------------------------------------------------------------------------
   procedure To_Hex(X : in Long_Long_Float; Digit_List : out String) is

      type Char_Array_T is array (Natural range <>) of Character;
      Hex_Array : constant Char_Array_T(0..15) :=
        ('0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F');

      Y : Long_Long_Float;
   begin
      Y := abs(X);
      for I in 1 .. Digit_List'Length loop
         Y := 16.0 * (Y-Long_Long_Float'Floor(Y));
         Digit_List(I) := Hex_Array(Integer(Long_Long_Float'Floor(Y)));
      end loop;
   end To_Hex;


   --------------------------------------------------------------------------------------
   --  Binary PI digit computation by the Bailey Borwein Plouffe formula:
   --
   --  Pi = SUM ( 1/(16 ** K) * ( 4/(8*K + 1) - 2/(8*K + 4) - 1/(8*K + 5) - 1/(8*K + 6))),
   --        Where K in [0, inf]
   --
   --------------------------------------------------------------------------------------
   procedure Bailey_Borwein_Plouffe(Digit : Natural) is

      S1 : constant Long_Long_Float := Evaluate_Serie(1.0, Digit-1);
      S2 : constant Long_Long_Float := Evaluate_Serie(4.0, Digit-1);
      S3 : constant Long_Long_Float := Evaluate_Serie(5.0, Digit-1);
      S4 : constant Long_Long_Float := Evaluate_Serie(6.0, Digit-1);

      Buffer_Length : constant Positive := 10;

      Position      : Long_Long_Float   := F_Mod((4.0 * S1) - (2.0 * S2) - S3 - S4, 1.0);
      Hex_Buffer    : String(1..Buffer_Length);

   begin
      if Position < 0.0 then
         Position := Position + 1.0;
      end if;

      To_Hex(Position, Hex_Buffer);
      Put_Line(Hex_Buffer);
   end Bailey_Borwein_Plouffe;


   --------------------------------------------------------------------------------------
   --
   --  PI digit computation variant by the Simon Plouffe
   --
   --------------------------------------------------------------------------------------
   function Plouffe(Digit : Natural) return Long_Long_Integer Is
      Av     : Long_Long_Integer;
      A      : Long_Long_Integer := 3;
      N      : Long_Long_Integer := Long_Long_Integer(Digit);
      NM     : Long_Long_Integer;
      DN     : Long_Long_Integer;
      K1     : Long_Long_Integer;
      K2     : Long_Long_Integer;
      T      : Long_Long_Integer;
      V      : Long_Long_Integer;
      V_Max  : Long_Long_Integer;
      S      : Long_Long_Integer;
      Sum    : Long_Long_Float := 0.0;

   begin
      N := Long_Long_Integer((Float(N +20) * Log(10.0,10.0))/Log(2.0,10.0));
      loop
         exit when A > (2*N);
         V_Max := Long_Long_Integer(Float'Floor(Log(Float(2*N),10.0)/Log(Float(A),10.0)));
         Av := 1;
         for I in 0 .. (V_Max - 1) loop
            Av := Av * A;
         end loop;

         S  := 0;
         NM := 1;
         DN := 1;
         V  := 0;
         K1 := 1;
         K2 := 1;

         for K in 1 .. N loop
            T := K;
            if K1 >= A then
               Inner_Loop: loop
                   T := Long_Long_Integer(Long_Long_Float(T)/Long_Long_Float(A));
                   V := V - 1;
                   exit Inner_Loop When T mod A /= 0;
               end loop Inner_Loop;
               K1 := 0;
            end if;
               K1 := K1 + 1;
               NM := (NM * T) mod Av;
               T := 2*K - 1;

               if K2 >= A then
                  if K2 = A then
                     Inner_Loop_1: loop
                       T := Long_Long_Integer(Float'Floor(Float(T/A)));
                       V := V + 1;
                       exit Inner_Loop_1 When T mod A /= 0;
                     end loop Inner_Loop_1;
                  end if;
                  K2 := K2 - A;
               end if;

               DN := (DN * T) mod Av;
               K2 := K2 + 2;

               if V > 0 then
                  T := Long_Long_Integer(Get_Inverse(Large_Positive(DN), Large_Positive(Av)));
                  T := (T * NM) mod Av;
                  T := (T * K) mod Av;

                  for I in V .. V_Max - 1 loop
                    T := (T * A) mod Av;
                  end loop;
                  S := S + T;
                  if S >= Av then
                     S := S - AV;
                  end if;
               end if;
            end loop;

           T := Long_Long_Integer(Slow_Reduce_Large_Exponent_Modulus(10, Positive(Digit-1), Large_Positive(Av)));
           S := (S*T) mod Av;
           Sum := F_Mod(Sum + Long_Long_Float(S) / Long_Long_Float(Av), 1.0);
           A:= Long_Long_Integer(Get_Prime_Slow(Large_Positive(A)));
      end loop;

      return Long_Long_Integer(Long_Long_Float'Floor(Sum *1.0E9));
   end Plouffe;
end PI_Calculation;
