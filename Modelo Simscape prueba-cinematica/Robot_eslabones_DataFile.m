% Simscape(TM) Multibody(TM) version: 7.4

% This is a model data file derived from a Simscape Multibody Import XML file using the smimport function.
% The data in this file sets the block parameter values in an imported Simscape Multibody model.
% For more information on this file, see the smimport function help page in the Simscape Multibody documentation.
% You can modify numerical values, but avoid any other changes to this file.
% Do not add code to this file. Do not edit the physical units shown in comments.

%%%VariableName:smiData


%============= RigidTransform =============%

%Initialize the RigidTransform structure array by filling in null values.
smiData.RigidTransform(11).translation = [0.0 0.0 0.0];
smiData.RigidTransform(11).angle = 0.0;
smiData.RigidTransform(11).axis = [0.0 0.0 0.0];
smiData.RigidTransform(11).ID = '';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(1).translation = [0 0 20];  % mm
smiData.RigidTransform(1).angle = 0;  % rad
smiData.RigidTransform(1).axis = [0 0 0];
smiData.RigidTransform(1).ID = 'B[pin2:1:-:Eslabon:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(2).translation = [-107.5 0 10];  % mm
smiData.RigidTransform(2).angle = 1.570796326794897;  % rad
smiData.RigidTransform(2).axis = [0 0 1];
smiData.RigidTransform(2).ID = 'F[pin2:1:-:Eslabon:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(3).translation = [107.5 0 10];  % mm
smiData.RigidTransform(3).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(3).axis = [0.70710678118654735 0.70710678118654768 0];
smiData.RigidTransform(3).ID = 'B[Eslabon:1:-:pin2:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(4).translation = [0 0 20];  % mm
smiData.RigidTransform(4).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(4).axis = [0 1 0];
smiData.RigidTransform(4).ID = 'F[Eslabon:1:-:pin2:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(5).translation = [107.5 0 10];  % mm
smiData.RigidTransform(5).angle = 1.5707963267948966;  % rad
smiData.RigidTransform(5).axis = [0 0 1];
smiData.RigidTransform(5).ID = 'B[Eslabon:2:-:pin1:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(6).translation = [0 0 10];  % mm
smiData.RigidTransform(6).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(6).axis = [0 0 1];
smiData.RigidTransform(6).ID = 'F[Eslabon:2:-:pin1:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(7).translation = [0 0 -10];  % mm
smiData.RigidTransform(7).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(7).axis = [0 0 1];
smiData.RigidTransform(7).ID = 'B[pin2:2:-:Eslabon:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(8).translation = [-107.5 0 -10];  % mm
smiData.RigidTransform(8).angle = 1.5707963267948968;  % rad
smiData.RigidTransform(8).axis = [-0 -0 -1];
smiData.RigidTransform(8).ID = 'F[pin2:2:-:Eslabon:2]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(9).translation = [0 0 0];  % mm
smiData.RigidTransform(9).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(9).axis = [1 0 0];
smiData.RigidTransform(9).ID = 'B[Soporte:1:-:pin2:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(10).translation = [0 0 0];  % mm
smiData.RigidTransform(10).angle = 3.1415926535897931;  % rad
smiData.RigidTransform(10).axis = [1 0 0];
smiData.RigidTransform(10).ID = 'F[Soporte:1:-:pin2:1]';

%Translation Method - Cartesian
%Rotation Method - Arbitrary Axis
smiData.RigidTransform(11).translation = [-8.1967211893347454 0.20628421448182249 -4.4999999999999982];  % mm
smiData.RigidTransform(11).angle = 0;  % rad
smiData.RigidTransform(11).axis = [0 0 0];
smiData.RigidTransform(11).ID = 'RootGround[Soporte:1]';


%============= Solid =============%
%Center of Mass (CoM) %Moments of Inertia (MoI) %Product of Inertia (PoI)

%Initialize the Solid structure array by filling in null values.
smiData.Solid(4).mass = 0.0;
smiData.Solid(4).CoM = [0.0 0.0 0.0];
smiData.Solid(4).MoI = [0.0 0.0 0.0];
smiData.Solid(4).PoI = [0.0 0.0 0.0];
smiData.Solid(4).color = [0.0 0.0 0.0];
smiData.Solid(4).opacity = 0.0;
smiData.Solid(4).ID = '';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(1).mass = 0.13793172641697282;  % kg
smiData.Solid(1).CoM = [0 -23.170056058546304 14.999999999999998];  % mm
smiData.Solid(1).MoI = [80.704880614624514 91.669251533384369 82.684947398327282];  % kg*mm^2
smiData.Solid(1).PoI = [0 0 0];  % kg*mm^2
smiData.Solid(1).color = [0.74901960784313726 0.74901960784313726 0.74901960784313726];
smiData.Solid(1).opacity = 1;
smiData.Solid(1).ID = 'Soporte.ipt_{A0606533-4230-6711-DFA6-8AAFFD99665E}';

%Inertia Type - Custom
%Visual Properties - Simple
smiData.Solid(2).mass = 0.0019006635554217435;  % kg
smiData.Solid(2).CoM = [1.1502474545203529e-09 0 10.000000000000002];  % mm
smiData.Solid(2).MoI = [0.077729220913435523 0.077729217200766684 0.028747534419419473];  % kg*mm^2
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
smiData.Solid(4).mass = 0.001335176877775605;  % kg
smiData.Solid(4).CoM = [1.0685923679086109e-09 0 2.6911764705882364];  % mm
smiData.Solid(4).MoI = [0.036475096332861688 0.036475094022846787 0.017886656618205293];  % kg*mm^2
smiData.Solid(4).PoI = [0 -8.1327113344396134e-13 0];  % kg*mm^2
smiData.Solid(4).color = [0.74901960784313726 0.74901960784313726 0.74901960784313726];
smiData.Solid(4).opacity = 1;
smiData.Solid(4).ID = 'pin1.ipt_{5FC51874-4119-D6DE-752A-91B152906F5F}';


%============= Joint =============%
%X Revolute Primitive (Rx) %Y Revolute Primitive (Ry) %Z Revolute Primitive (Rz)
%X Prismatic Primitive (Px) %Y Prismatic Primitive (Py) %Z Prismatic Primitive (Pz) %Spherical Primitive (S)
%Constant Velocity Primitive (CV) %Lead Screw Primitive (LS)
%Position Target (Pos)

%Initialize the RevoluteJoint structure array by filling in null values.
smiData.RevoluteJoint(2).Rz.Pos = 0.0;
smiData.RevoluteJoint(2).ID = '';

smiData.RevoluteJoint(1).Rz.Pos = 9.086900000000373;  % deg
smiData.RevoluteJoint(1).ID = '[pin2:1:-:Eslabon:1]';

smiData.RevoluteJoint(2).Rz.Pos = -18.169999999999931;  % deg
smiData.RevoluteJoint(2).ID = '[pin2:2:-:Eslabon:2]';

