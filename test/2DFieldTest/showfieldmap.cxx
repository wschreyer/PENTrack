int showfieldmap(){
	mytree->Draw("sqrt(Bx*Bx+By*By+Bz*Bz) - abs(cos(sqrt(x*x+y*y)*pi) + cos(z*pi)):sqrt(x*x+y*y):z","","");
}
