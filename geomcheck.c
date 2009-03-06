#include "main.h"

int GeomCheck(long double *y,long double r,long double vr,long double z,long double vz,long double phi,long double vphi, long double t)
{
	int AlreadyReflected = 0, InVolume=0;
	long double rNorm,zNorm;
	
	
	// Absorberhit
	if((r>=absrmin) && (r<=absrmax))
	{
		if((z>=abszmin) && (z<=abszmax))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			AlreadyReflected = Absorber(y,y[1],y[2],y[3],y[4],y[5],y[6],t);
		}
	}
	
	
	// round outer bottom corner
	if( (r>=RoundBottomCornerCenterr-wanddicke)&&(z>=RoundBottomCornerCenterz-wanddicke)&&(sqrtl(powl(r-RoundBottomCornerCenterr,2)+powl(z-RoundBottomCornerCenterz,2))<RoundBottomCornerradius) )
	{
		InVolume++;	
		if(t==0) return 1;    // if initial point is already invalid return 1
			rNorm = (r-RoundBottomCornerCenterr), zNorm = z-RoundBottomCornerCenterz; // normalenvektor der Kreisoberfläche
			AlreadyReflected =ReflexGeneral(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado, rNorm,zNorm,0.0);
	}
	
	// Storage volume
	if((r>=StorVolrmin-wanddicke)&&(r<=StorVolrmax+wanddicke)&&(z>=StorVolzmin-wanddicke)&&(z<=StorVolzmax))
	{
		InVolume++;	
		// reflection at BOTTOM, but not if channel or round corner is there, but do reflect at channel blockage
		if((y[4]<0)&&(!stopall)&&(z<=StorVolzmin)&&(r<StorVolrmax)&&( ( (r<RoundBottomCornerCenterr)&&(r<FillChannelrmin) )||(fabsl(fmod((double) (phi/conv), 90.0))<(FillChannelBlockageAngle/2)) )     )
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			//AlreadyReflected =ReflexBottom(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
			AlreadyReflected =ReflexBottom(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealDLC,FPimDLC); 
		}		
		// reflection at INNER wall
		if((y[2]<0)&&(!stopall)&&(r<=StorVolrmin))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			//AlreadyReflected =ReflexIn(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
			AlreadyReflected =ReflexIn(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealDLC,FPimDLC);
		}		
		// reflection at OUTER wall
		if((y[2]>0)&&(!stopall)&&(r>=StorVolrmax)&&( (z>=FillChannelzmax)||(fabsl(fmod((double) (phi/conv),90.0))<(FillChannelBlockageAngle/2)) )  )
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			//AlreadyReflected =ReflexOut(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
			AlreadyReflected =ReflexOut(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealDLC,FPimDLC);
		}		
		
		// configuration with a bigger detector at 1.15 m height and no focusing coils
		if(detz<StorVolzmax)
		{
				// DETECTOR hit 
		if((z>=detz)&&(z<=detz+0.01)&&(r>=detrmin)&&(r<=detrmax))
		{
			if (y[4]>0)     //Treffer auf den Detektor von unten
				{
					if(t==0) return 1;    // if initial point is already invalid return 1
					AlreadyReflected =ReflexTop(y,y[1],y[2],y[3],y[4],y[5],y[6],t,0.5,FPrealCsI,FPimCsI);  // reflection on CsI, rough -> probability of diffuse reflection 0.5
					if(stopall)
					{
						kennz=KENNZAHL_DETECTOR_BOTTOM;        // 6
						if(vend>0) gammaend= acosl(y[4]/vend) /conv;
						else gammaend=0;
						alphaend= atan2l(y[6]*r,y[2])/conv;			
						return 0;						
					}
				
				}
				
				if (y[4]>0)     //Treffer auf den Detektor von oben
				{
					if(t==0) return 1;    // if initial point is already invalid return 1
					AlreadyReflected =ReflexBottom(y,y[1],y[2],y[3],y[4],y[5],y[6],t,0.5,FPrealCsI,FPimCsI);  // reflection on CsI, rough -> probability of diffuse reflection 0.5
					if(stopall)
					{
						kennz=KENNZAHL_DETECTOR_BOTTOM;        // 6
						if(vend>0) gammaend= acosl(y[4]/vend) /conv;
						else gammaend=0;
						alphaend= atan2l(y[6]*r,y[2])/conv;			
						return 0;						
					}
				
				}
			
		}
		// escape to TOP
		if((y[4]>0)&&(z>=StorVolzmax))
			{                          //Entkommen nach oben kennz=KENNZAHL_ESCAPE_TOP
				if(t==0) return 1;    // if initial point is already invalid return 1
				AlreadyReflected =ReflexTop(y,y[1],y[2],y[3],y[4],y[5],y[6],t,0.5,FPrealNocado,FPimNocado);  // reflection on top -> has to be isolator -> PE is assumed
				if(stopall)
				{		
					if(vend>0) gammaend= acosl(y[4]/vend) /conv;
					else gammaend=0;
					alphaend= atan2l(y[6]*r,y[2])/conv;
					kennz=KENNZAHL_ESCAPE_TOP;    // 4
					stopall=1;
					return 0;
				}
			}
		}
	}
	
	// filling channel
	if((r>=FillChannelrmin-wanddicke)&&(r<=FillChannelrmax+wanddicke)&&(z>=FillChannelzmin-wanddicke)&&(z<=FillChannelzmax+wanddicke))
	{
		// no BOTTOM
		InVolume++;	
		// reflection at INNER wall, but not if still in storage volume or at round corner
		if((y[2]<0)&&(!stopall)&&(r<=FillChannelrmin)&&(z<StorVolzmin)&&(z>Bufferzmax)&&(z<RoundBottomCornerCenterz))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			AlreadyReflected = ReflexIn(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
		}		
		// reflection at OUTER wall
		if((y[2]>0)&&(!stopall)&&(r>=FillChannelrmax))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			AlreadyReflected = ReflexOut(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
		}				
		// reflection at CONE
		if((z>=(FillConezbot+(FillConerbot-r)/(FillConertop-FillConerbot)*(FillConeztop-FillConezbot)))   && 	(r>=(FillConerbot-(z-FillConezbot)/(FillConeztop-FillConezbot)*fabsl(FillConertop-FillConerbot))) && (r>=FillConertop)	&& (z>=FillConezbot))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			rNorm = -(FillConeztop-FillConezbot), zNorm = FillConertop-FillConerbot; // normalenvektor des Konus
			AlreadyReflected =ReflexGeneral(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado, rNorm,zNorm,0.0);
		}		
		
				
	}
	
	// buffer volume
	if((r>=Bufferrmin-wanddicke)&&(r<=Bufferrmax+wanddicke)&&(z>=Bufferzmin-wanddicke)&&(z<=Bufferzmax+wanddicke))
	{
		InVolume++;	
		// reflection at BOTTOM, filling tube at center, and UCN detector
		if(z<=Bufferzmin)
		{
			// UCN not in entrance tube nor UCN det
			if((y[4]<0)&&(!stopall)&& ((r>UCNentrancermax)||(!slit))&&( (sqrtl( powl(r*cosl(phi)-UCNdetr*cosl(UCNdetphi),2) + powl(r*sinl(phi)-UCNdetr*sinl(UCNdetphi),2) )>UCNdetradius)||(!DetOpen) )  )
			{
				if(t==0) return 1;    // if initial point is already invalid return 1
				AlreadyReflected = ReflexBottom(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
			}
			// UCN detector hit => 100% efficiency is assumed
			if((y[4]<0)&&(DetOpen)&& ( sqrtl( powl(r*cosl(phi)-UCNdetr*cosl(UCNdetphi),2) + powl(r*sinl(phi)-UCNdetr*sinl(UCNdetphi),2) )  <UCNdetradius) )
			{
				kennz=12;  // UCN detector hit
				stopall=1;
				return 0;
			}
			// UCN in entrance tube => also loss
			if((y[4]<0)&&(r<=UCNentrancermax)&&(slit))
			{
				kennz=5;        // particle in entrance hole                                                                          
				stopall=1;
				return 0;
			}
		}		
		// reflection at OUTER wall
		if((y[2]>0)&&(!stopall)&&(r>=Bufferrmax))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			AlreadyReflected = ReflexOut(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
		}		
		// reflection at TOP wall 
		if((y[4]>0)&&(!stopall)&&(z>=Bufferzmax)&&( (r<=FillChannelrmin)||(fabsl(fmod((double) (phi/conv), 90.0))<(FillChannelBlockageAngle/2)) ))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			AlreadyReflected =ReflexTop(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
		}		
		
	}
	
	// detector volume
	if((r>=DetVolrmin-wanddicke)&&(r<=DetVolrmax+wanddicke)&&(z>=DetVolzmin-wanddicke)&&(z<=DetVolzmax+wanddicke))
	{
		
		InVolume++;	
		// reflection at BOTTOM: inside of storvolrmin
		if((y[4]<0)&&(!stopall)&&(z<=DetVolzmin)&&(r<StorVolrmin))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			AlreadyReflected =ReflexBottom(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
		}		
		
		// reflection at INNER wall
		if((y[2]<0)&&(!stopall)&&(r<=DetVolrmin))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			AlreadyReflected =ReflexIn(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
		}		
		// reflection at OUTER wall
		if((y[2]>0)&&(!stopall)&&(r>=DetVolrmax))
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			AlreadyReflected =ReflexOut(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado);
		}				
		// DETECTOR hit 
		if((z>=detz)&&(z<=detz+0.01)&&(r>=detrmin)&&(r<=detrmax))
		{
			if (y[4]>0)     //Treffer auf den Detektor von unten
				{
					if(t==0) return 1;    // if initial point is already invalid return 1
					AlreadyReflected =ReflexTop(y,y[1],y[2],y[3],y[4],y[5],y[6],t,0.5,FPrealCsI,FPimCsI);  // reflection on CsI, rough -> probability of diffuse reflection 0.5
					if(stopall)
					{
						kennz=KENNZAHL_DETECTOR_BOTTOM;        // 6
						if(vend>0) gammaend= acosl(y[4]/vend) /conv;
						else gammaend=0;
						alphaend= atan2l(y[6]*r,y[2])/conv;			
						return 0;						
					}
				
				}
			
		}
		// escape to TOP
		if((y[4]>0)&&(z>=DetVolzmax))
			{                          //Entkommen nach oben kennz=KENNZAHL_ESCAPE_TOP
				if(t==0) return 1;    // if initial point is already invalid return 1
				AlreadyReflected =ReflexTop(y,y[1],y[2],y[3],y[4],y[5],y[6],t,0.5,FPrealPE,FPimPE);  // reflection on top -> has to be isolator -> PE is assumed
				if(stopall)
				{		
					if(vend>0) gammaend= acosl(y[4]/vend) /conv;
					else gammaend=0;
					alphaend= atan2l(y[6]*r,y[2])/conv;
					kennz=KENNZAHL_ESCAPE_TOP;    // 4
					stopall=1;
					return 0;
				}
			}
		// reflection at CONE
		if((z>=(DetConezbot+(DetConerbot-r)/(DetConertop-DetConerbot)*(DetConeztop-DetConezbot)))   && 	r>=(DetConerbot-(z-DetConezbot)/(DetConeztop-DetConezbot)*fabsl(DetConertop-DetConerbot))	)
		{
			if(t==0) return 1;    // if initial point is already invalid return 1
			rNorm = -(DetConeztop-DetConezbot), zNorm = DetConertop-DetConerbot; // normalenvektor des Konus			
			AlreadyReflected =ReflexGeneral(y,y[1],y[2],y[3],y[4],y[5],y[6],t,DiffProb,FPrealNocado,FPimNocado, rNorm,zNorm,0.0);
		}
			
	}
	
	if(InVolume==0)
	{
		return 1;
	}
	
		
	return 0;
}

