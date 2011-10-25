-------------------------------------------------------------------------------
--                                                                           --
--                                   PI                                      --
--                                                                           --
--                                 PI.adb                                    --
--                                                                           --
--                                  BODY                                     --
--                                                                           --
--                   Copyright (C) 1996 Ulrik Hørlyk Hjort                   --
--                                                                           --
--  PI is free software;  you can  redistribute it                           --
--  and/or modify it under terms of the  GNU General Public License          --
--  as published  by the Free Software  Foundation;  either version 2,       --
--  or (at your option) any later version.                                   --
--  PI is distributed in the hope that it will be                            --
--  useful, but WITHOUT ANY WARRANTY;  without even the  implied warranty    --
--  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                  --
--  See the GNU General Public License for  more details.                    --
--  You should have  received  a copy of the GNU General                     --
--  Public License  distributed with Yolk.  If not, write  to  the  Free     --
--  Software Foundation,  51  Franklin  Street,  Fifth  Floor, Boston,       --
--  MA 02110 - 1301, USA.                                                    --
--                                                                           --
-------------------------------------------------------------------------------
with PI_Calculation; use PI_Calculation;
with Ada.Text_IO; use Ada.Text_IO;

procedure PI is

   Ret_Val_Plouffe : Long_Long_Integer;

begin
   Monte_Carlo_Method(1.0E-5);
   Gregory_Leibniz_Series(1.0E-6);
   Wallis_Product(1.0E-6);
   New_Line;
   Put_Line("Digits calculated by Plouffe's formular");
   Put("Digit pos: 3 = ");
   Bailey_Borwein_Plouffe(3);
   Put_Line("Digits calculated by Plouffe's formular");
   Put("Digit pos: 5 = ");
   Bailey_Borwein_Plouffe(5);
   Put_Line("Digits calculated by Plouffe's formular");
   Put("Digit pos: 30 = ");
   Bailey_Borwein_Plouffe(30);
   Put_Line("Digits calculated by Plouffe's formular");
   Put("Digit pos: 151 = ");
   Bailey_Borwein_Plouffe(151);
   New_Line;
   New_Line;
   Put_Line("Digits calculated by Plouffe's formular");
   Ret_Val_Plouffe := Plouffe(3);
   Put_Line("Digit pos: 3 = " & " " & Long_Long_Integer'Image(Ret_Val_Plouffe));
   Ret_Val_Plouffe := Plouffe(5);
   Put_Line("Digit pos: 5 = " & " " & Long_Long_Integer'Image(Ret_Val_Plouffe));
   Ret_Val_Plouffe := Plouffe(30);
   Put_Line("Digit pos: 30 = " & " " & Long_Long_Integer'Image(Ret_Val_Plouffe));
   Ret_Val_Plouffe := Plouffe(151);
   Put_Line("Digit pos: 151 = " & " " & Long_Long_Integer'Image(Ret_Val_Plouffe));
end PI;
