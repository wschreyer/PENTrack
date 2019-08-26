void showfieldmap(){

	// Drawing inherited from 2DFieldTest
	// mytree->Draw("sqrt(Bx*Bx+By*By+Bz*Bz) - abs(cos(sqrt(x*x+y*y)*pi) + cos(z*pi)):sqrt(x*x+y*y):z","","");
	
	// Basic Drawing of Bz magnitude through a cut in z
	mytree->Draw("Bz:x","","");
}
