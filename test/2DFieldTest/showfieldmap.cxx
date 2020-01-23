void showfieldmap(){
	TCanvas *c1 = new TCanvas("c1", "c1");
	mytree->Draw("sqrt(Bx**2 + By**2 + Bz**2) - abs(cos(sqrt(x*x+y*y)*pi) + cos(z*pi)):sqrt(x*x+y*y):z","","");
	TCanvas *c2 = new TCanvas("c2", "c2");
	mytree->Draw("abs(V/10000) - abs(cos(sqrt(x*x+y*y)*pi) + cos(z*pi)):sqrt(x*x+y*y):z","","");
	TCanvas *c3 = new TCanvas("c3", "c3");
	mytree->Draw("sqrt(Ex**2+Ey**2+Ez**2)/10000 - pi*sqrt(sin(pi*sqrt(x**2+y**2))**2 + sin(pi*z)**2):z:sqrt(x**2+y**2)");
}
