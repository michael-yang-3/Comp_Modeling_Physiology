% Michael Yang | my2699
% Due Mar 07, 2021
% BMEN E6003
%% Homework 2C Problem 1 - 2 element Hill Model

% to run your code simply hit "Run" or use the shortcut "command(Mac)/
% control(Windows)+enter"

clear;clc

% Create  your  own  two-element* Hill  model  in  MATLAB,  assuming  length
% L(t)  is prescribed  as  an  input.  Confirm  that  your  implementation  
% gives  expected  force responses  for  the  development  of  isometric  
% tension,  step  release,  and  constant velocity release.
% *Include Xse and Xce in your code.

% 2 Elements = Xse + Xce


% Case 1 - ISOMETRIC TENSION

% Setup L and t
L = [repmat(1,1000,1)];
t = [linspace(0, 5, 1000)'];

% Hill.p to get value for Lse
[P,H,Lse,Lce] = hill(L,t);

% Setup force
p = [repmat(0, 1000, 1)];

%constants
a = 399 * 0.098;
b = 0.331;
p_0 = a / 0.235;
alpha = (p_0 / Lse(1));

% Loop
for i=1:999
    Lse(i) = Lse(1) + (p(i) / alpha);
    Lce(i) = L(i) - Lse(i);
    dt=t(i+1)-t(i);
    dL=L(i+1)-L(i);
    dp = alpha * [(dL/dt) + ((p_0 - p(i)) / (p(i) + a))*b]*dt;
    p(i+1)=p(i)+(dp);
end
% 
figure(1)
hold on; subplot(2,2,1); plot(t, p,'k');axis([0 5 0 200]); title('Force'); xlabel('Time (sec)'); ylabel('Force (mN/mm2)');
hold on; subplot(2,2,2); plot(t, L, 'r');axis([0 5 0 2]); title('Length'); xlabel('Time (sec)'); ylabel('Length (mm)');
hold on; subplot(2,2,3); plot(t, Lse, 'b'); axis([0 5 0 1]); title('Lse'); xlabel('Time (sec)'); ylabel('Length (mm)'); 
hold on; subplot(2,2,4); plot(t, Lce, 'r'); axis([0 5 0 1]); title('Lce'); xlabel('Time (sec)'); ylabel('Length (mm)');


% Case 2 - Step Change Release

% Setup L and t
L2 = [repmat(1,400,1);repmat(0.9,600,1)];
t = [linspace(0, 5, 1000)'];

% Hill.p to get values for Lse
[P,H,Lse2,Lce2] = hill(L2,t);

% Setup force
p = [repmat(0, 1000, 1)];

%constants
a = 399 * 0.098;
b = 0.331;
p_0 = a / 0.235;
alpha = (p_0 / Lse2(1));

% Loop
for i=1:999
    Lse2(i) = Lse2(1) + (p(i) / alpha);
    Lce2(i) = L2(i) - Lse2(i);
    dt=t(i+1)-t(i);
    dL2=L2(i+1)-L2(i);
    dp2 = alpha * [(dL2/dt) + ((p_0 - p(i)) / (p(i) + a))*b]*dt;
    p(i+1)=p(i)+(dp2);
end
% 
figure(2)
hold on; subplot(2,2,1); plot(t, p,'k');axis([0 5 0 200]); title('Force'); xlabel('Time (sec)'); ylabel('Force (mN/mm2)');
hold on; subplot(2,2,2); plot(t, L2, 'r');axis([0 5 0 2]); title('Length'); xlabel('Time (sec)'); ylabel('Length (mm)');
hold on; subplot(2,2,3); plot(t, Lse2, 'b'); axis([0 5 0 1]); title('Lse'); xlabel('Time (sec)'); ylabel('Length (mm)'); 
hold on; subplot(2,2,4); plot(t, Lce2, 'r'); axis([0 5 0 1]); title('Lce'); xlabel('Time (sec)'); ylabel('Length (mm)');



% Case 3 - Constant Velocity Release

% Setup L and t
L3 = [repmat(1,300,1);linspace(1,0.85,700)'];
t = [linspace(0, 5, 1000)'];

% Hill.p to find Lse value
[P,H,Lse3,Lce3] = hill(L3,t);

% Setup force
p = [repmat(0, 1000, 1)];

%constants
a = 399 * 0.098;
b = 0.331;
p_0 = a / 0.235;
alpha = (p_0 / Lse3(1));

% Loop
for i=1:999
    Lse3(i) = Lse3(1) + (p(i) / alpha);
    Lce3(i) = L3(i) - Lse3(i);
    dt=t(i+1)-t(i);
    dL3=L3(i+1)-L3(i);
    dp3 = alpha * [(dL3/dt) + ((p_0 - p(i)) / (p(i) + a))*b]*dt;
    p(i+1)=p(i)+(dp3);
end
% 
figure(3)
hold on; subplot(2,2,1); plot(t, p,'k');axis([0 5 0 170]); title('Force'); xlabel('Time (sec)'); ylabel('Force (mN/mm2)');
hold on; subplot(2,2,2); plot(t, L3, 'r');axis([0 5 0 2]); title('Length'); xlabel('Time (sec)'); ylabel('Length (mm)');
hold on; subplot(2,2,3); plot(t, Lse3, 'b'); axis([0 5 0 1]); title('Lse'); xlabel('Time (sec)'); ylabel('Length (mm)'); 
hold on; subplot(2,2,4); plot(t, Lce3, 'r'); axis([0 5 0 1]); title('Lce'); xlabel('Time (sec)'); ylabel('Length (mm)');

