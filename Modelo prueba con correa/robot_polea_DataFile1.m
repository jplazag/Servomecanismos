% Simscape(TM) Multibody(TM) version: 7.4

% This is a model data file derived from a Simscape Multibody Import XML file using the smimport function.
% The data in this file sets the block parameter values in an imported Simscape Multibody model.
% For more information on this file, see the smimport function help page in the Simscape Multibody documentation.
% You can modify numerical values, but avoid any other changes to this file.
% Do not add code to this file. Do not edit the physical units shown in comments.

%%%VariableName:smiData


%============= RigidTransform =============%

%Initialize the RigidTransform structure array by filling in null values.
smiData.RigidTransform(21).translation = [0.0 0.0 0.0];
smiData.RigidTransform(21).angle = 0.0;
smiData.RigidTransform(21).axis = [0.0 0.0 0.0];
smiData.RigidTransform(21).ID = '';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(1).translation = [0 0 -20];  % mm
smiData.RigidTransform(1).angle = 1.5707963267948968;  % rad
smiData.RigidTransform(1).axis = [-0 -0 -1];
smiData.RigidTransform(1).ID = 'B[Soporte:1:-:Eslabon:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(2).translation = [-107.5 0 0];  % mm
smiData.RigidTransform(2).angle = 1.5707963267948968;  % rad
smiData.RigidTransform(2).axis = [-0 -0 -1];
smiData.RigidTransform(2).ID = 'F[Soporte:1:-:Eslabon:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(3).translation = [107.5 0 0];  % mm
smiData.RigidTransform(3).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(3).axis = [0.70710678118654746 0.70710678118654757 0];
smiData.RigidTransform(3).ID = 'B[Eslabon:1:-:pin2:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(4).translation = [0 0 0];  % mm
smiData.RigidTransform(4).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(4).axis = [0 1 0];
smiData.RigidTransform(4).ID = 'F[Eslabon:1:-:pin2:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(5).translation = [0 0 20];  % mm
smiData.RigidTransform(5).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(5).axis = [0 0 1];
smiData.RigidTransform(5).ID = 'B[pin2:2:-:Eslabon:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(6).translation = [-107.5 0 10];  % mm
smiData.RigidTransform(6).angle = 1.5707963267948968;  % rad
smiData.RigidTransform(6).axis = [0 0 1];
smiData.RigidTransform(6).ID = 'F[pin2:2:-:Eslabon:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(7).translation = [0 0 -30];  % mm
smiData.RigidTransform(7).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(7).axis = [-0.70710678118654746 0.70710678118654757 0];
smiData.RigidTransform(7).ID = 'B[Soporte:1:-:pin2:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(8).translation = [0 0 0];  % mm
smiData.RigidTransform(8).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(8).axis = [0 1 0];
smiData.RigidTransform(8).ID = 'F[Soporte:1:-:pin2:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(9).translation = [0 0 0];  % mm
smiData.RigidTransform(9).angle = 0;  % rad
smiData.RigidTransform(9).axis = [0 0 0];
smiData.RigidTransform(9).ID = 'B[V-Belt transmission:1:-:]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(10).translation = [0 0 0];  % mm
smiData.RigidTransform(10).angle = 0;  % rad
smiData.RigidTransform(10).axis = [0 0 0];
smiData.RigidTransform(10).ID = 'F[V-Belt transmission:1:-:]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(11).translation = [0 0 0];  % mm
smiData.RigidTransform(11).angle = 0;  % rad
smiData.RigidTransform(11).axis = [0 0 0];
smiData.RigidTransform(11).ID = 'B[V-Belt transmission:1:-:pin2:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(12).translation = [-0.72267823673903031 12.636611833557733 44.999999999999986];  % mm
smiData.RigidTransform(12).angle = 1.5707963267948968;  % rad
smiData.RigidTransform(12).axis = [0 0 1];
smiData.RigidTransform(12).ID = 'F[V-Belt transmission:1:-:pin2:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(13).translation = [0 0 0];  % mm
smiData.RigidTransform(13).angle = 0;  % rad
smiData.RigidTransform(13).axis = [0 0 0];
smiData.RigidTransform(13).ID = 'B[V-Belt transmission:1:-:pin2:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(14).translation = [165.48598588305538 117.59254543750906 34.999999999999986];  % mm
smiData.RigidTransform(14).angle = 2.9172672103921631;  % rad
smiData.RigidTransform(14).axis = [-0 -0 -1];
smiData.RigidTransform(14).ID = 'F[V-Belt transmission:1:-:pin2:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(15).translation = [0 0.75000000000000178 0];  % mm
smiData.RigidTransform(15).angle = 0;  % rad
smiData.RigidTransform(15).axis = [0 0 0];
smiData.RigidTransform(15).ID = 'B[Soporte:1:-:]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(16).translation = [-12.636611833557733 0.027321763260973464 -14.999999999999982];  % mm
smiData.RigidTransform(16).angle = 0;  % rad
smiData.RigidTransform(16).axis = [0 0 0];
smiData.RigidTransform(16).ID = 'F[Soporte:1:-:]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(17).translation = [-12.636611833557733 0.027321763260973464 -11.062999999999983];  % mm
smiData.RigidTransform(17).angle = 0;  % rad
smiData.RigidTransform(17).axis = [0 0 0];
smiData.RigidTransform(17).ID = 'AssemblyGround[V-Belt transmission:1:V-Belt:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(18).translation = [-12.636611833557733 -0.72267823673902654 -16.062999999999988];  % mm
smiData.RigidTransform(18).angle = 0;  % rad
smiData.RigidTransform(18).axis = [0 0 0];
smiData.RigidTransform(18).ID = 'AssemblyGround[V-Belt transmission:1:Grooved Pulley1:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(19).translation = [187.49794804898229 77.834039748550978 -16.062999999999988];  % mm
smiData.RigidTransform(19).angle = 0;  % rad
smiData.RigidTransform(19).axis = [0 0 0];
smiData.RigidTransform(19).ID = 'AssemblyGround[V-Belt transmission:1:Grooved Pulley2:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(20).translation = [-12.636611833557733 -0.72267823673902831 -14.999999999999982];  % mm
smiData.RigidTransform(20).angle = 0;  % rad
smiData.RigidTransform(20).axis = [0 0 0];
smiData.RigidTransform(20).ID = 'RootGround[Soporte:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(21).translation = [0 0 0];  % mm
smiData.RigidTransform(21).angle = 0;  % rad
smiData.RigidTransform(21).axis = [0 0 0];
smiData.RigidTransform(21).ID = 'RootGround[V-Belt transmission:1]';


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
smiData.Solid(1).color = [0.74901960784313726 0.74901960784313726 0.74901960784313726];
smiData.Solid(1).opacity = 1;
smiData.Solid(1).ID = 'Soporte.ipt_{A0606533-4230-6711-DFA6-8AAFFD99665E}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(2).mass = 0.0038013271108434869;  % kg
smiData.Solid(2).CoM = [1.1502474545203529e-09 0 20.000000000000004];  % mm
smiData.Solid(2).MoI = [0.5355911529112195 0.5355911454858816 0.057495068838838946];  % kg*mm^2
smiData.Solid(2).PoI = [0 0 0];  % kg*mm^2
smiData.Solid(2).color = [0.74901960784313726 0.74901960784313726 0.74901960784313726];
smiData.Solid(2).opacity = 1;
smiData.Solid(2).ID = 'pin2.ipt_{04912817-4149-7EDD-C109-1CA4B5F7B0D8}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(3).mass = 0.096665707058937247;  % kg
smiData.Solid(3).CoM = [1.133743443233888e-11 0 5];  % mm
smiData.Solid(3).MoI = [13.514477521573633 479.51588805103745 491.4192704549622];  % kg*mm^2
smiData.Solid(3).PoI = [0 0 0];  % kg*mm^2
smiData.Solid(3).color = [0.74901960784313726 0.74901960784313726 0.74901960784313726];
smiData.Solid(3).opacity = 1;
smiData.Solid(3).ID = 'Eslabon.ipt_{AE7C2FE2-4E44-6EDB-EB86-C1AC78484A19}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(4).mass = 0.0055159831902384234;  % kg
smiData.Solid(4).CoM = [100.06727994126999 38.528358992645025 0];  % mm
smiData.Solid(4).MoI = [4.1333284943191728 23.987269876190041 28.104613258689142];  % kg*mm^2
smiData.Solid(4).PoI = [0 0 -9.2124331544267015];  % kg*mm^2
smiData.Solid(4).color = [0.25098039215686274 0.25098039215686274 0.25098039215686274];
smiData.Solid(4).opacity = 1;
smiData.Solid(4).ID = 'V-Belt_{707F1982-4EF7-6323-1265-0F8F978203D0}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(5).mass = 0.0230594366289412;  % kg
smiData.Solid(5).CoM = [0 3.227267409805407e-09 4.9999999999999956];  % mm
smiData.Solid(5).MoI = [0.89494950725709965 0.89494976690746941 1.2855003401623579];  % kg*mm^2
smiData.Solid(5).PoI = [0 0 0];  % kg*mm^2
smiData.Solid(5).color = [0.8784313725490196 0.87450980392156863 0.85882352941176465];
smiData.Solid(5).opacity = 1;
smiData.Solid(5).ID = 'Grooved Pulley1_{DF90438A-4F86-11A9-A893-14A928022253}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(6).mass = 0.0230594366289412;  % kg
smiData.Solid(6).CoM = [0 3.227267409805407e-09 4.9999999999999956];  % mm
smiData.Solid(6).MoI = [0.89494950725709965 0.89494976690746941 1.2855003401623579];  % kg*mm^2
smiData.Solid(6).PoI = [0 0 0];  % kg*mm^2
smiData.Solid(6).color = [0.8784313725490196 0.87450980392156863 0.85882352941176465];
smiData.Solid(6).opacity = 1;
smiData.Solid(6).ID = 'Grooved Pulley2_{581F8C68-49DF-4AA8-30A3-F0B910F44E67}';


%============= Joint =============%
%X Revolute Primitive (Rx) %Y Revolute Primitive (Ry) %Z Revolute Primitive (Rz)
%X Prismatic Primitive (Px) %Y Prismatic Primitive (Py) %Z Prismatic Primitive (Pz) %Spherical Primitive (S)
%Constant Velocity Primitive (CV) %Lead Screw Primitive (LS)
%Position Target (Pos)

%Initialize the PlanarJoint structure array by filling in null values.
smiData.PlanarJoint(1).Rz.Pos = 0.0;
smiData.PlanarJoint(1).Px.Pos = 0.0;
smiData.PlanarJoint(1).Py.Pos = 0.0;
smiData.PlanarJoint(1).ID = '';

smiData.PlanarJoint(1).Rz.Pos = 0;  % deg
smiData.PlanarJoint(1).Px.Pos = 0;  % mm
smiData.PlanarJoint(1).Py.Pos = 0;  % mm
smiData.PlanarJoint(1).ID = '[Soporte:1:-:]';


%Initialize the RevoluteJoint structure array by filling in null values.
smiData.RevoluteJoint(3).Rz.Pos = 0.0;
smiData.RevoluteJoint(3).ID = '';

smiData.RevoluteJoint(1).Rz.Pos = 21.430975568420411;  % deg
smiData.RevoluteJoint(1).ID = '[Soporte:1:-:Eslabon:1]';

%This joint has been chosen as a cut joint. Simscape Multibody treats cut joints as algebraic constraints to solve closed kinematic loops. The imported model does not use the state target data for this joint.
smiData.RevoluteJoint(2).Rz.Pos = 124.28387670104628;  % deg
smiData.RevoluteJoint(2).ID = '[Eslabon:1:-:pin2:2]';

smiData.RevoluteJoint(3).Rz.Pos = 179.99999999999997;  % deg
smiData.RevoluteJoint(3).ID = '[Soporte:1:-:pin2:1]';

