import graph;

size(8cm);

real xmin = 0;
real xmax = 5;
real ymin = -2.5;
real ymax = 2.5;

real w0 = 2.5;
real w = 1;

real fg(real x) { return 2.5*x^2/(5^2); } 
real fl(real x) { return -2.5*x^2/(5^2); }
real f0(real x) { return 0;}

real Fg(real x) { return -2*2.5*x/(5^2);} 
real Fl(real x) { return 2*2.5*x/(5^2);}


pair Bg = (xmax, fg(xmax));
pair Bl = (xmax, fl(xmax));
pair B0 = (xmax, f0(xmax));

xaxis(xmin, xmax);
yaxis(ymin, ymax);

draw(graph(fg, xmin, xmax), red+linewidth(1));
draw(graph(fl, xmin, xmax), green+linewidth(1));
draw(graph(f0, xmin, xmax), blue+linewidth(1));

draw((xmin+1, 1)--(xmin+0.2, 1),  arrow=Arrow(), red);
draw((xmin+0.2,-1)--(xmin+1, -1), arrow=Arrow(), green);


draw((xmin,0)--(xmax,0),Arrow);
draw((0,ymin)--(0,ymax),Arrow);



label("$\rho$", (xmax,0), E);
label("$U_{\rm eff}$", (0,ymax), N);


label("$\omega_0>\omega$", Bg, NE, red);
label("$\omega_0<\omega$", Bl, NE, green);
label("$\omega_0=\omega$", B0, NW, blue);