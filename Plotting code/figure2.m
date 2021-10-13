% Produce the panels in Figure 2 of our manuscript.

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/), and requires the export_fig package
% (freely available at https://github.com/altmany/export_fig) to export the
% figure panels to pdf files.

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/gen_tost_serial_indep.mat','f_gen_indep','f_tost_indep','f_serial_indep')
load('../Results/gen_tost_serial_mech.mat','f_gen_mech','f_tost_mech','f_serial_mech')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Pre-format figures

for k = 1:3
figsetup(k)
end

% Generation time

figure(1); hold on;
p1 = plot(f_gen_indep,'color',c1,'linewidth',2);
p2 = plot(f_gen_mech,'color',c2,'linewidth',2);
xlim([0,15])
ylim([0,0.25])
xlabel('Generation time (days)')
ylabel('Density')
l = legend([p1,p2],{'Independent transmission/symptoms','Mechanistic'});
l.FontSize = 15;
l.Position = [0.2200    0.8110    0.5800    0.0810];

% Time from onset of symptoms to transmission (TOST)

figure(2); hold on;
plot(f_tost_indep,'color',c1,'linewidth',2)
plot(f_tost_mech,'color',c2,'linewidth',2)
xlim([-10,10])
ylim([0,0.2])
xlabel({'Time from onset of symptoms';'to transmission (days)'})
ylabel('Density')
l = legend('Independent transmission/symptoms','Mechanistic');
l.FontSize = 15;
l.Position = [0.2200    0.8110    0.5800    0.0810];

% Serial interval

figure(3); hold on;
plot(f_serial_indep,'color',c1,'linewidth',2)
plot(f_serial_mech,'color',c2,'linewidth',2)
xlim([-5,20])
xlabel('Serial interval (days)')
ylabel('Density')
l = legend('Independent transmission/symptoms','Mechanistic');
l.FontSize = 15;
l.Position = [0.2200    0.8110    0.5800    0.0810];

% Mean and standard deviation of serial interval

t_indep = chebfun('t',f_serial_indep.domain);
t_mech = chebfun('t',f_serial_mech.domain);

m_ser_indep = sum(t_indep*f_serial_indep);
s_ser_indep = sqrt(sum(t_indep*t_indep*f_serial_indep)-m_ser_indep^2);
m_ser_mech = sum(t_mech*f_serial_mech);
s_ser_mech = sqrt(sum(t_mech*t_mech*f_serial_mech)-m_ser_mech^2);

% Post-format figures

for k = 1:3
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/figure2/figure2A.pdf -painters -nocrop -transparent
figure(2); export_fig Figures/figure2/figure2B.pdf -nocrop -transparent
figure(3); export_fig Figures/figure2/figure2C.pdf -nocrop -transparent

rmpath('../Data')
rmpath('../Results')