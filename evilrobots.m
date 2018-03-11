alpha=[-pi/2,0,pi/2,-pi/2,pi/2,0]; % Twist angles (in radians)
l=[0,0.5,0,0,0,0]; % Link lengths a (in metres), called l to distinguish from the a component of DH matrix
d=[0,0.25,0,1,0,0.5]; % Offset distances (in metres)
T06=[-sqrt(2)/2,0,sqrt(2)/2,1;0,-1,0,1;sqrt(2)/2,0,sqrt(2)/2,0;0,0,0,1]; % Composite Denavit-Hartenburg matrix
n=T06(1:3,1);
s=T06(1:3,2);
a=T06(1:3,3);
p=T06(1:3,4);
parm=p-d(6)*a;
pxarm=parm(1);
pyarm=parm(2);
pzarm=parm(3);
ax=a(1);
ay=a(2);
az=a(3);
sx=s(1);
sy=s(2);
sz=s(3);
nx=n(1);
ny=n(2);
nz=n(3);
% Joint angles
t1l=atan2(pyarm,pxarm)-atan2(d(2),-sqrt(pxarm^2+pyarm^2-d(2)^2)); % theta1 for left arm/shoulder configuration
t1r=atan2(pyarm,pxarm)-atan2(d(2),sqrt(pxarm^2+pyarm^2-d(2)^2)); % theta1 for right arm/shoulder configuration
theta1=[t1l,t1r];
A=cos(theta1)*pxarm+sin(theta1)*pyarm;
B=(A.^2+pzarm^2+l(2)^2-d(4)^2)/(2*l(2));
t2u=atan2(A,pzarm)-atan2(B,-sqrt(A.^2+pzarm^2-B.^2)); % theta2 for elbow up configuration
t2d=atan2(A,pzarm)-atan2(B,sqrt(A.^2+pzarm^2-B.^2)); % theta2 for elbow down configuration
theta2=[t2u;t2d];
theta1=[theta1;theta1]; % Duplicating theta 1 and A to make them 2x2 matrices so that matrix dimensions agree in calculations
A=[A;A];
theta3=atan2(A-l(2)*cos(theta2),pzarm+l(2)*sin(theta2))-theta2;
theta4=atan2(-sin(theta1)*ax+cos(theta1)*ay,cos(theta2+theta3).*(cos(theta1)*ax+sin(theta1)*ay)-sin(theta2+theta3)*az);
theta5=atan2(sqrt((cos(theta1).*cos(theta2+theta3)*ax+sin(theta1).*cos(theta2+theta3)*ay-sin(theta2+theta3)*az).^2+(-sin(theta1)*ax+cos(theta1)*ay).^2),sin(theta2+theta3).*(cos(theta1)*ax+sin(theta1)*ay)+cos(theta2+theta3)*az);
theta6=atan2(sin(theta2+theta3).*(cos(theta1)*sx+sin(theta1)*sy)+cos(theta2+theta3)*sz,-(sin(theta2+theta3).*(cos(theta1)*nx+sin(theta1)*ny)+cos(theta2+theta3)*nz));
% Forward kinematics
% A matrices
theta=zeros(1,4,6); % Pre-allocating memory
theta(:,:,1)=[theta1(1),theta1(2),theta1(3),theta1(4)]; % Creating 3D 1x4x6 theta matrix to allow index to range from 1 to 4 rather than 1 to 2 twice and save two extra nested for loops below
theta(:,:,2)=[theta2(1),theta2(2),theta2(3),theta2(4)];
theta(:,:,3)=[theta3(1),theta3(2),theta3(3),theta3(4)];
theta(:,:,4)=[theta4(1),theta4(2),theta4(3),theta4(4)];
theta(:,:,5)=[theta5(1),theta5(2),theta5(3),theta5(4)];
theta(:,:,6)=[theta6(1),theta6(2),theta6(3),theta6(4)];
A=zeros(4,4,4,6);
for u=linspace(1,4,4) % Calculating A matrices for each of the 4 sets of thetas = 24 matrices, 2 of the A01s will be duplicates as 2 duplicate theta1 values need to be fed in
    for v=linspace(1,6,6) % Calculating 6 A matrices for each theta
        A(:,:,u,v)=[cos(theta(:,u,v)),-cos(alpha(v))*sin(theta(:,u,v)),sin(alpha(v))*sin(theta(:,u,v)),l(v)*cos(theta(:,u,v));sin(theta(:,u,v)),cos(alpha(v))*cos(theta(:,u,v)),-sin(alpha(v))*cos(theta(:,u,v)),l(v)*sin(theta(:,u,v));0,sin(alpha(v)),cos(alpha(v)),d(v);0,0,0,1];
    end
end
A01=A(:,:,:,1);
A12=A(:,:,:,2);
A23=A(:,:,:,3);
A34=A(:,:,:,4);
A45=A(:,:,:,5);
A56=A(:,:,:,6);
T1=zeros(4,4,4);
T2=zeros(4,4,4);
T=zeros(4,4,4);
for t=linspace(1,4,4)
    T1(:,:,t)=A01(:,:,t)*A12(:,:,t)*A23(:,:,t); % T1 matrices = A03
    T2(:,:,t)=A34(:,:,t)*A45(:,:,t)*A56(:,:,t); % T2 matrices = A36
    T(:,:,t)=T1(:,:,t)*T2(:,:,t); % Composite DH matrix T = T1*T2
end
R03=T1(1:3,1:3,:); % Rotation matrices R03
R36=T2(1:3,1:3,:); % Rotation matrices R36
R06=T(1:3,1:3,:); % Rotation matrices R06
% User output selection
outputinverse=0;
outputforward=0;
outputtrans=0;
outputDH=0;
outputrot=0;
skipselection=0;
kinematics=input('Would you like the results for the forward or inverse kinematics? Type forward, inverse or both.\nAlternatively type skip to skip result selection and output all results at once.\n','s'); % Get user input
if strcmp(kinematics,'both') % Compare user input string to expected input and output only the desired results as otherwise output can be long and messy due to having 24 (really 22) A matrices
    outputinverse=1;
    outputforward=1;
elseif strcmp(kinematics,'inverse')
    outputinverse=1;
elseif strcmp(kinematics,'forward')
    outputforward=1;
elseif strcmp(kinematics,'skip')
    outputinverse=1;
    outputforward=1;
    outputtrans=1;
    outputDH=1;
    outputrot=1;
    skipselection=1;
else
    fprintf('Input error, please type "inverse" to request i) theta solutions, "forward" to request ii) matrix solutions, "both" to request both sets of solutions at once or "skip" (without "") to skip selection and output all solutions.\n')
    return
end
if skipselection==0
    if outputforward==1
        trans=input('Would you like the homogeneous transformation matrices A01-A56 for each set of joint angles? Type yes or no.\n','s');
        DH=input('Would you like the Denavit Hartenburg transformation matrix T=T1*T2 for each set of joint angles? Type yes or no.\n','s');
        rot=input('Would you like the rotation matrices R03, R36 and R06 for each set of joint angles? Type yes or no.\n','s');
        if strcmp(trans,'yes')
            outputtrans=1;
        elseif strcmp(trans,'no')
        else
            fprintf('Input error, please type "yes" or "no" (without "") to request the homogeneous transformation matrices A01-A56.\n')
            return % End program if input doesn't match either yes or no
        end
        if strcmp(DH,'yes')
            outputDH=1;
        elseif strcmp(DH,'no')
        else
            fprintf('Input error, please type "yes" or "no" to request the Denavit Hartenburg transformation matrix T=T1*T2.\n')
            return
        end
        if strcmp(rot,'yes')
            outputrot=1;
        elseif strcmp(rot,'no')
        else
            fprintf('Input error, please type "yes" or "no" to request the rotation matrices R03, R36 and R06.\n')
            return
        end
    end
end
% Output
theta1=[t1l,t1r]; % Changing theta1 back to its real format for output
if outputinverse==1
    fprintf('The inverse kinematic solutions are shown below with the joint angles in degrees, with theta 1 displayed as a matrix in the format:\n')
    fprintf('left arm       right arm\n')
    fprintf('Thetas 2-6 are displayed as matrices in the format:\n')
    fprintf('left arm elbow up      right arm elbow up\nleft arm elbow down     right arm elbow down\n')
    fprintf('The matrix element each joint angle appears in corresponds to its elbow and arm configuration,\n')
    fprintf('i.e. the upper-right element of each matrix gives the joint angle for the elbow up, right arm configuration, etc.\n')
    fprintf('Theta 1\n')
    disp(theta1*180/pi)
    fprintf('Theta 2\n')
    disp(theta2*180/pi)
    fprintf('Theta 3\n')
    disp(theta3*180/pi)
    fprintf('Theta 4\n')
    disp(theta4*180/pi)
    fprintf('Theta 5\n')
    disp(theta5*180/pi)
    fprintf('Theta 6\n')
    disp(theta6*180/pi)
end
if outputforward==1
    fprintf('The forward kinematic solutions are:\n')
    for w=linspace(1,4,4) % Outputs the following for each set of joint angles
        if outputtrans==1 || outputDH==1 || outputrot==1
            if w==1
                fprintf('For the left arm, elbow up configuration,\n')
            elseif w==2
                fprintf('For the right arm, elbow up configuration,\n')
            elseif w==3
                fprintf('For the left arm, elbow down configuration,\n')
            else
                fprintf('For the right arm, elbow down configuration,\n')
            end
        end
        if outputtrans==1
            fprintf('The homogeneous transformation matrices A01-A56 are:\n')
            fprintf('A01\n') % The two A01 right matrices will be identical to each other, as will the two A01 left matrices but they are kept to easily output the solutions for the whole configuration as w ranges from 1 to 4
            disp(A01(:,:,w))
            fprintf('A12\n')
            disp(A12(:,:,w))
            fprintf('A23\n')
            disp(A23(:,:,w))
            fprintf('A34\n')
            disp(A34(:,:,w))
            fprintf('A45\n')
            disp(A45(:,:,w))
            fprintf('A56\n')
            disp(A56(:,:,w))
        end
        if outputDH==1
            fprintf('The Denavit-Hartenburg transformation matrix T is:\n')
            disp(T(:,:,w))
        end
        if outputrot==1
            fprintf('The Rotation matrix R03 is:\n')
            disp(R03(:,:,w))
            fprintf('The Rotation matrix R36 is:\n')
            disp(R36(:,:,w))
            fprintf('The Rotation matrix R06 is:\n')
            disp(R06(:,:,w))
        end
    end
tol=1e-15;
if ismembertol(T(:,:,1),T(:,:,2),tol) & ismembertol(T(:,:,2),T(:,:,3),tol) & ismembertol(T(:,:,3),T(:,:,4),tol) %#ok<AND2> % Checking whether matrices are identical to the chosen tolerance to verify that forward kinematics have worked identically for each set of joint angles
    fprintf('The Denavit-Hartenburg transformation matrices T are identical for whichever arm and elbow configuration is used to calculate T to a %g tolerance (due to floating point representation).\n',tol)
    if ismembertol(T(:,:,1),T06,tol) % Checking whether matrices are identical to the original T06 matrix to verify that forward kinematics have regained original T06 matrix
        disp(T06)
        fprintf('They are also identical to a %g tolerance to the T06 matrix shown above which was given in the question. This verifies that the calculated inverse kinematic joint angles are correct.\n',tol)
    end
else
    fprintf('Something has gone horribly wrong.\n')
end
end