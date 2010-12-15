#ifndef RACETRACK_H_
#define RACETRACK_H_

struct TRacetrack{
	long double I_rt;
	TRacetrack(long double I): I_rt(I){};
	virtual ~TRacetrack(){};
	virtual void BFeld(long double x, long double y, long double z, long double B[4][4]) = 0;
};

struct TFiniteWire: public TRacetrack{
	long double SW1x,SW1y,SW1z, SW2x,SW2y,SW2z;
	TFiniteWire(long double SW1xx, long double SW1yy, long double SW1zz, long double SW2xx, long double SW2yy, long double SW2zz, long double I): TRacetrack(I), SW1x(SW1xx), SW1y(SW1yy), SW1z(SW1zz), SW2x(SW2xx), SW2y(SW2yy), SW2z(SW2zz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

struct TFiniteWireX: public TRacetrack{
	long double SW1x,SW2x,SWz;
	TFiniteWireX(long double SW1xx, long double SW2xx, long double SWzz, long double I): TRacetrack(I), SW1x(SW1xx), SW2x(SW2xx), SWz(SWzz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

struct TFiniteWireY: public TRacetrack{
	long double SW1y,SW2y,SWz;
	TFiniteWireY(long double SW1yy, long double SW2yy, long double SWzz, long double I): TRacetrack(I), SW1y(SW1yy), SW2y(SW2yy), SWz(SWzz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

struct TFiniteWireZ: public TRacetrack{
	long double SWx,SWy,SW1z,SW2z;
	TFiniteWireZ(long double SWxx, long double SWyy, long double SW1zz, long double SW2zz, long double I): TRacetrack(I), SWx(SWxx), SWy(SWyy), SW1z(SW1zz), SW2z(SW2zz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

struct TFiniteWireZCenter: public TRacetrack{
	long double SW1z, SW2z;
	TFiniteWireZCenter(long double SW1zz, long double SW2zz, long double I): TRacetrack(I), SW1z(SW1zz), SW2z(SW2zz){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

struct TFullRacetrack: public TRacetrack{
	long double SW1z, SW2z, SWr;
	TFullRacetrack(long double SW1zz, long double SW2zz, long double SWrr, long double I): TRacetrack(I), SW1z(SW1zz), SW2z(SW2zz), SWr(SWrr){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

struct TInfiniteWireZ: public TRacetrack{
	long double lx, ly;
	TInfiniteWireZ(long double lxx, long double lyy, long double I): TRacetrack(I), lx(lxx), ly(lyy){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

struct TInfiniteWireZCenter: public TRacetrack{
	TInfiniteWireZCenter(long double I): TRacetrack(I){};
	void BFeld(long double x, long double y, long double z, long double B[4][4]);
};

#endif /*RACETRACK_H_*/
