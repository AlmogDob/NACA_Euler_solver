clc; clear; close all


fig1 = figure ("Name","J vals mat",'Position',[100 300 900 500]);
J_vasl_mat = readmatrix("checks\J_vals_mat.txt");
m = mesh(J_vasl_mat);
m.FaceColor='interp';
title("J vals mat");
colorbar

fig2 = figure ("Name","dxi dx mat",'Position',[200 300 900 500]);
dxi_dx_mat = readmatrix("checks\dxi_dx.txt");
m = mesh(dxi_dx_mat);
m.FaceColor='interp';
title("dxi dx mat");
colorbar

fig3 = figure ("Name","dxi dy mat",'Position',[300 300 900 500]);
dxi_dy_mat = readmatrix("checks\dxi_dy.txt");
m = mesh(dxi_dy_mat);
m.FaceColor='interp';
title("dxi dy mat");
colorbar

fig4 = figure ("Name","deta dx mat",'Position',[400 300 900 500]);
deta_dx_mat = readmatrix("checks\deta_dx.txt");
m = mesh(deta_dx_mat);
m.FaceColor='interp';
title("deta dx mat");
colorbar

fig5 = figure ("Name","deta dy mat",'Position',[500 300 900 500]);
deta_dy_mat = readmatrix("checks\deta_dy.txt");
m = mesh(deta_dy_mat);
m.FaceColor='interp';
title("deta dy mat");
colorbar

fig6 = figure ("Name","dx dxi mat",'Position',[100 50 900 500]);
dx_dxi_mat = readmatrix("checks\dx_dxi.txt");
m = mesh(dx_dxi_mat);
m.FaceColor='interp';
title("dx dxi mat");
colorbar

fig7 = figure ("Name","dy dxi mat",'Position',[200 50 900 500]);
dy_dxi_mat = readmatrix("checks\dy_dxi.txt");
m = mesh(dy_dxi_mat);
m.FaceColor='interp';
title("dy dxi mat");
colorbar

fig8 = figure ("Name","dx deta mat",'Position',[300 50 900 500]);
dx_deta_mat = readmatrix("checks\dx_deta.txt");
m = mesh(dx_deta_mat);
m.FaceColor='interp';
title("dx deta mat");
colorbar

fig9 = figure ("Name","dy deta mat",'Position',[400 50 900 500]);
dy_deta_mat = readmatrix("checks\dy_deta.txt");
m = mesh(dy_deta_mat);
m.FaceColor='interp';
title("dy deta mat");
colorbar

fig10 = figure ("Name","experimant",'Position',[400 50 900 500]);

m = mesh(J_vasl_mat.*dy_deta_mat);
m.FaceColor='interp';
title("experimant");
colorbar


