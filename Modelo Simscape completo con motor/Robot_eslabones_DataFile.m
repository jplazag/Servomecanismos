% Simscape(TM) Multibody(TM) version: 7.4

% This is a model data file derived from a Simscape Multibody Import XML file using the smimport function.
% The data in this file sets the block parameter values in an imported Simscape Multibody model.
% For more information on this file, see the smimport function help page in the Simscape Multibody documentation.
% You can modify numerical values, but avoid any other changes to this file.
% Do not add code to this file. Do not edit the physical units shown in comments.

%%%VariableName:smiData


%============= RigidTransform =============%

%Initialize the RigidTransform structure array by filling in null values.
smiData.RigidTransform(19).translation = [0.0 0.0 0.0];
smiData.RigidTransform(19).angle = 0.0;
smiData.RigidTransform(19).axis = [0.0 0.0 0.0];
smiData.RigidTransform(19).ID = '';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(1).translation = [13.215473463303709 0 0];  % mm
smiData.RigidTransform(1).angle = 1.9235588160682087;  % rad
smiData.RigidTransform(1).axis = [0.50673185397138654 -0.69745656233305509 -0.50673185397138654];
smiData.RigidTransform(1).ID = 'B[Polea:1:-:Correa:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(2).translation = [0 0 5.25];  % mm
smiData.RigidTransform(2).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(2).axis = [0.47843583569762133 0.87812251486926285 0];
smiData.RigidTransform(2).ID = 'F[Polea:1:-:Correa:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(3).translation = [13.215473463303709 0 0];  % mm
smiData.RigidTransform(3).angle = 1.9235588160682087;  % rad
smiData.RigidTransform(3).axis = [0.50673185397138654 -0.69745656233305497 -0.50673185397138654];
smiData.RigidTransform(3).ID = 'B[Polea:2:-:Correa:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(4).translation = [215 0 5.25];  % mm
smiData.RigidTransform(4).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(4).axis = [0.47843583569762133 0.87812251486926285 0];
smiData.RigidTransform(4).ID = 'F[Polea:2:-:Correa:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(5).translation = [0 0 20];  % mm
smiData.RigidTransform(5).angle = 0;  % rad
smiData.RigidTransform(5).axis = [0 0 0];
smiData.RigidTransform(5).ID = 'B[pin2:1:-:Eslabon:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(6).translation = [-107.5 0 10];  % mm
smiData.RigidTransform(6).angle = 1.5707963267948961;  % rad
smiData.RigidTransform(6).axis = [0 0 1];
smiData.RigidTransform(6).ID = 'F[pin2:1:-:Eslabon:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(7).translation = [107.5 0 10];  % mm
smiData.RigidTransform(7).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(7).axis = [0.70710678118654768 0.70710678118654735 0];
smiData.RigidTransform(7).ID = 'B[Eslabon:1:-:pin2:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(8).translation = [0 0 20];  % mm
smiData.RigidTransform(8).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(8).axis = [0 1 0];
smiData.RigidTransform(8).ID = 'F[Eslabon:1:-:pin2:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(9).translation = [107.5 0 10];  % mm
smiData.RigidTransform(9).angle = 1.5707963267948968;  % rad
smiData.RigidTransform(9).axis = [0 0 1];
smiData.RigidTransform(9).ID = 'B[Eslabon:2:-:pin1:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(10).translation = [0 0 10];  % mm
smiData.RigidTransform(10).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(10).axis = [0 0 1];
smiData.RigidTransform(10).ID = 'F[Eslabon:2:-:pin1:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(11).translation = [0 0 -10];  % mm
smiData.RigidTransform(11).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(11).axis = [0 0 1];
smiData.RigidTransform(11).ID = 'B[pin2:2:-:Eslabon:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(12).translation = [-107.5 0 -10];  % mm
smiData.RigidTransform(12).angle = 1.5707963267948968;  % rad
smiData.RigidTransform(12).axis = [-0 -0 -1];
smiData.RigidTransform(12).ID = 'F[pin2:2:-:Eslabon:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(13).translation = [0 0 -30];  % mm
smiData.RigidTransform(13).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(13).axis = [1 0 0];
smiData.RigidTransform(13).ID = 'B[Soporte:1:-:Polea:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(14).translation = [13.398500000000002 0 0];  % mm
smiData.RigidTransform(14).angle = 2.8342931083834495;  % rad
smiData.RigidTransform(14).axis = [-0.69857538902427807 -0.15487043520038671 0.69857538902427807];
smiData.RigidTransform(14).ID = 'F[Soporte:1:-:Polea:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(15).translation = [0 0 0];  % mm
smiData.RigidTransform(15).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(15).axis = [0 0 1];
smiData.RigidTransform(15).ID = 'B[pin2:2:-:Polea:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(16).translation = [13.398500000000004 0 0];  % mm
smiData.RigidTransform(16).angle = 3.0059423100619784;  % rad
smiData.RigidTransform(16).axis = [-0.70547345838991893 0.067929367837001778 -0.70547345838991893];
smiData.RigidTransform(16).ID = 'F[pin2:2:-:Polea:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(17).translation = [0 0 -30];  % mm
smiData.RigidTransform(17).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(17).axis = [1 0 0];
smiData.RigidTransform(17).ID = 'B[Soporte:1:-:pin2:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(18).translation = [0 0 0];  % mm
smiData.RigidTransform(18).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(18).axis = [1 0 0];
smiData.RigidTransform(18).ID = 'F[Soporte:1:-:pin2:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(19).translation = [-8.1967211893347454 0.20628421448182249 -4.4999999999999982];  % mm
smiData.RigidTransform(19).angle = 0;  % rad
smiData.RigidTransform(19).axis = [0 0 0];
smiData.RigidTransform(19).ID = 'RootGround[Soporte:1]';


%============= Solid =============%
%Center of Mass (CoM) %Moments of Inertia (MoI) %Product of Inertia (PoI)

%Initialize the Solid structure array by filling in null values.
smiData.Solid(6).mass = 0.0;
smiData.Solid(6).CoM = [0.0 0.0 0.0];
smiData.Solid(6).MoI = [0.0 0.0 0.0];
smiData.Solid(6).PoI = [0.0 0.0 0.0];
smiData.Solid(6).color = [0.0 0.0 0.0];
smiData.Solid(6).opacity = 0.0;
smiData.Solid(6).ID = '';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(1).mass = 0.036128781874027724;  % kg
smiData.Solid(1).CoM = [3.02634936765929e-11 -1.5068780025989061 -14.999999999999998];  % mm
smiData.Solid(1).MoI = [9.5492217973447548 10.101728743557036 11.823047801529123];  % kg*mm^2
smiData.Solid(1).PoI = [0 0 -1.6390685220807174e-12];  % kg*mm^2
smiData.Solid(1).color = [0.49803921568627452 0.49803921568627452 0.49803921568627452];
smiData.Solid(1).opacity = 1;
smiData.Solid(1).ID = 'Soporte.ipt_{A0606533-4230-6711-DFA6-8AAFFD99665E}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(2).mass = 0.00076969020012946656;  % kg
smiData.Solid(2).CoM = [7.3197617197742704e-10 0 10.000000000000002];  % mm
smiData.Solid(2).MoI = [0.02801351639442359 0.028013515785577325 0.0047143521713698582];  % kg*mm^2
smiData.Solid(2).PoI = [0 0 0];  % kg*mm^2
smiData.Solid(2).color = [0.74901960784313726 0.74901960784313726 0.74901960784313726];
smiData.Solid(2).opacity = 1;
smiData.Solid(2).ID = 'pin2.ipt_{04912817-4149-7EDD-C109-1CA4B5F7B0D8}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(3).mass = 0.13936026959027722;  % kg
smiData.Solid(3).CoM = [2.8829542716112931e-12 0 5];  % mm
smiData.Solid(3).MoI = [19.288684416976785 701.96517486875041 718.93118812588921];  % kg*mm^2
smiData.Solid(3).PoI = [0 0 0];  % kg*mm^2
smiData.Solid(3).color = [0.45098039215686275 0.27843137254901962 0.039215686274509803];
smiData.Solid(3).opacity = 0.55000001192092896;
smiData.Solid(3).ID = 'Eslabon.ipt_{AE7C2FE2-4E44-6EDB-EB86-C1AC78484A19}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(4).mass = 0.00054650622203070121;  % kg
smiData.Solid(4).CoM = [6.7911355920512801e-10 0 2.6017964071856294];  % mm
smiData.Solid(4).MoI = [0.013068052249618094 0.013068051867740458 0.0029569134476357204];  % kg*mm^2
smiData.Solid(4).PoI = [0 -2.1669352273955777e-13 0];  % kg*mm^2
smiData.Solid(4).color = [0.74901960784313726 0.74901960784313726 0.74901960784313726];
smiData.Solid(4).opacity = 1;
smiData.Solid(4).ID = 'pin1.ipt_{5FC51874-4119-D6DE-752A-91B152906F5F}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(5).mass = 0.0023018718543368296;  % kg
smiData.Solid(5).CoM = [8.7181716846870376 -0.017148367814004715 -0.017047726591491372];  % mm
smiData.Solid(5).MoI = [0.087048467397501342 0.11360577943497062 0.11360462531759784];  % kg*mm^2
smiData.Solid(5).PoI = [6.79462283886543e-07 0.0003421394665713602 0.00034416356490770441];  % kg*mm^2
smiData.Solid(5).color = [0.77254901960784317 0.76078431372549016 0.58039215686274515];
smiData.Solid(5).opacity = 1;
smiData.Solid(5).ID = 'Polea.ipt_{B69E30E1-4B47-5049-0BC4-838C7D69ACAE}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(6).mass = 0.0085858323348192017;  % kg
smiData.Solid(6).CoM = [107.50000000000281 2.3089365776817544e-13 8.0812780218861421e-14];  % mm
smiData.Solid(6).MoI = [0.74889869983233226 42.109541768109139 42.700675798789177];  % kg*mm^2
smiData.Solid(6).PoI = [0 1.8560353254315682e-13 -0.0003566679525017343];  % kg*mm^2
smiData.Solid(6).color = [0.25098039215686274 0.25098039215686274 0.25098039215686274];
smiData.Solid(6).opacity = 1;
smiData.Solid(6).ID = 'Correa.ipt_{D6036E33-45E7-9EA4-83B1-EAA46AF1B844}';


%============= Joint =============%
%X Revolute Primitive (Rx) %Y Revolute Primitive (Ry) %Z Revolute Primitive (Rz)
%X Prismatic Primitive (Px) %Y Prismatic Primitive (Py) %Z Prismatic Primitive (Pz) %Spherical Primitive (S)
%Constant Velocity Primitive (CV) %Lead Screw Primitive (LS)
%Position Target (Pos)

%Initialize the RevoluteJoint structure array by filling in null values.
smiData.RevoluteJoint(5).Rz.Pos = 0.0;
smiData.RevoluteJoint(5).ID = '';

smiData.RevoluteJoint(1).Rz.Pos = 41.484695995388741;  % deg
smiData.RevoluteJoint(1).ID = '[Polea:1:-:Correa:1]';

smiData.RevoluteJoint(2).Rz.Pos = -137.51719812094561;  % deg
smiData.RevoluteJoint(2).ID = '[Polea:2:-:Correa:1]';

smiData.RevoluteJoint(3).Rz.Pos = 152.68189144719778;  % deg
smiData.RevoluteJoint(3).ID = '[pin2:1:-:Eslabon:1]';

smiData.RevoluteJoint(4).Rz.Pos = -48.71328267892477;  % deg
smiData.RevoluteJoint(4).ID = '[pin2:2:-:Eslabon:2]';

%This joint has been chosen as a cut joint. Simscape Multibody treats cut joints as algebraic constraints to solve closed kinematic loops. The imported model does not use the state target data for this joint.
smiData.RevoluteJoint(5).Rz.Pos = 14.316214436472135;  % deg
smiData.RevoluteJoint(5).ID = '[pin2:2:-:Polea:2]';

