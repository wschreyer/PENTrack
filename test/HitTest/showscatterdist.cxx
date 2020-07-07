int showscatterdist(){
	std::string phi1 = "atan2(v1x,v1z)";
	std::string phi2 = "atan2(v2x,v2z)";
	std::string theta1 = "acos((v1x*nx+v1y*ny+v1z*nz)/sqrt(v1x*v1x+v1y*v1y+v1z*v1z))";
	std::string theta2 = "acos((v2x*nx+v2y*ny+v2z*nz)/sqrt(v2x*v2x+v2y*v2y+v2z*v2z))";
	int total = neutronhit->GetEntries();
	int reflected = neutronhit->GetEntries((theta2 + " < TMath::Pi()/2").c_str());
	int transmitted = neutronhit->GetEntries((theta2 + " >= TMath::Pi()/2").c_str());
	int specularReflected = neutronhit->GetEntries(("(" + theta1 + " - TMath::Pi()/2 - " + theta2 + ")**2 + (" + phi1 + " - " + phi2 + ")**2 < 1e-6").c_str());
	int specularTransmitted = neutronhit->GetEntries((theta2 + " >= TMath::Pi()/2 && (" + phi1 + " - " + phi2 + ")**2 < 1e-16").c_str());
	std::cout << "Reflected: " << (double)reflected/total << " +/- " << std::sqrt(reflected)/total << " (Expected: 0.597797)" << std::endl;
	std::cout << "Transmitted: " << (double)transmitted/total << " +/- " << std::sqrt(transmitted)/total << " (Expected: 0.402203)" << std::endl;
	std::cout << "Specularly reflected: " << (double)specularReflected/reflected << " +/- " << std::sqrt(specularReflected)/reflected << " (Expected: 0.9)" << std::endl;
	std::cout << "Specularly transmitted: " << (double)specularTransmitted/transmitted << " +/- " << std::sqrt(specularTransmitted)/transmitted << " (Expected: 0.9)" << std::endl;
	return neutronhit->Draw((phi2 + ':' + theta2).c_str(),"","COLZ");
}
