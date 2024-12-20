

closeNoPrompt(matlab.desktop.editor.getAll);
close all
desktop = com.mathworks.mde.desk.MLDesktop.getInstance();
 desktop.closeGroup('Variables');
clear
clc