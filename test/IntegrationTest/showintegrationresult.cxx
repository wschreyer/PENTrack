void showintegrationresult(){
	// a = -g - mu*dBdz/m
	// z(t) = z_start + v_z_start*t + a/2*t^2
	// t(z=zend) = (-v_z_start +- sqrt(v_z_start^2 - 2*a*(z_start - z_end)))/a
	// x_hit = x_start + v_x_start*t(z=0)
	neutronend->Draw("xend - xstart - vxstart*(-vzstart - sqrt(vzstart*vzstart - 2*-15.5754762729*(zstart - zend)))/-15.5754762729");
}
