int showscatterdist(){
	return neutronhit->Draw("atan2(v2x,v2z):acos((v2x*nx+v2y*ny+v2z*nz)/sqrt(v2x*v2x+v2y*v2y+v2z*v2z))","","COLZ");
}
