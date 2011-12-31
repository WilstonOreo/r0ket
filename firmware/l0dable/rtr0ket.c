// This file is part of rtr0ket. rtr0ket is free software: you can
// redistribute it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, version 2.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License aint with
// this program; if not, write to the Free Software Foundation, Inc., 51
// Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// Copyright (2011), Wilston Oreo

#define L0DABLE

#include <sysinit.h>

#include "basic/basic.h"
#include "lcd/render.h"
#include "lcd/display.h"
#include "lcd/allfonts.h"

#ifdef L0DABLE
#include "usetable.h"
void main_rtr0ket();

void ram(void)
{
	main_rtr0ket();
}
#endif
#ifndef L0DABLE
#include "../l0dable/usetable.h"
#endif

#define XSIZE 96L
#define YSIZE 68L

typedef unsigned char u8;

const u8 logo[] = { 0x00,0x00,0xee,0x2a,
	0xee,0x8a,0xee,0x00,
	0xee,0x82,0x86,0x82,
	0xee,0x00,0x00,0x00 };
const u8 pesthorn [] = { 
	0x00,0x00,
	0x00,0x00,
	0x47,0xce,
	0x6f,0xee,
	0x79,0x2e,
	0x71,0x1e,
	0x30,0x1c,
	0x3e,0xfc,
	0x1f,0xf8,
	0x0a,0xa0,
	0x60,0x06,
	0x76,0xde,
	0x6c,0x66,
	0x00,0x00,
	0x00,0x00,
	0x00,0x00};

u8 nicknameBuf [XSIZE/8*8];

typedef int Veci3[3];
int posX = 65536;
int posY = 256;
int timer = 0;
int angle = 0;
int speed = 0;
int mode  = 0;

///////// Look up tables for required for integer sin() and integer cos()

const u8 sin_table[] = {
	0, 13, 25, 38, 50, 62, 75, 86, 98, 110, 121, 132,
	142, 152, 162, 172, 181, 189, 198, 205, 213, 219, 
	225, 231, 236, 241, 245, 248, 251, 253, 254, 255, 255 };

int intsin(int x)
{
	int a,b,sn;
	if (x < 0) x = -x -256; x &= 1023;
	if (x & 256)
	{
		int pos = 31 - ((x/8) & 31);

		a = sin_table[pos];
		b = sin_table[pos+1];
		sn = a + (b-a)*(7-(x & 7))/8; 
	} else
	{
		int pos = (x/8) & 31;
		a = sin_table[pos];
		b = sin_table[pos+1];
		sn = a + (b-a)*(x & 7)/8;
	}
	if (x & 512) sn = -sn;
	return sn;
}

int intcos(int x)
{
	return intsin(x+256);
}

/* by Mark Crowne 
 * http://www.azillionmonkeys.com/qed/sqroot.html
 * */

unsigned int intsqrt (unsigned int val) {
	unsigned int temp, g=0;

	if (val >= 0x40000000L) {
		g = 0x8000; 
		val -= 0x40000000L;
	}

#define INNER_ISQRT(s)                        \
	temp = (g << (s)) + (1 << ((s) * 2 - 2));   \
	if (val >= temp) {                          \
		g += 1 << ((s)-1);                        \
		val -= temp;                              \
	}

	INNER_ISQRT (15)
		INNER_ISQRT (14)
		INNER_ISQRT (13)
		INNER_ISQRT (12)
		INNER_ISQRT (11)
		INNER_ISQRT (10)
		INNER_ISQRT ( 9)
		INNER_ISQRT ( 8)
		INNER_ISQRT ( 7)
		INNER_ISQRT ( 6)
		INNER_ISQRT ( 5)
		INNER_ISQRT ( 4)
		INNER_ISQRT ( 3)
		INNER_ISQRT ( 2)
#undef INNER_ISQRT
		temp = g+g+1;
	if (val >= temp) g++;
	return g;
}

// Ray-sphere intersection
int sphereIntersection(Veci3 spherePos, Veci3 rayOrg, Veci3 rayDir, int R)
{
	int t = -1000; 
	Veci3 ocVec;
	Veci3 rayDiri;
	ocVec[0] = (rayOrg[0] - spherePos[0]) >> 8;
	ocVec[1] = (rayOrg[1] - spherePos[1]) >> 8;
	ocVec[2] = (rayOrg[2] - spherePos[2]) >> 8;
	rayDiri[0] = rayDir[0] >> 8;
	rayDiri[1] = rayDir[1] >> 8;
	rayDiri[2] = rayDir[2] >> 8;
	// The equation to solve has the form
	// At² + Bt + C
	int A = rayDiri[0]*rayDiri[0] + rayDiri[1]*rayDiri[1] + rayDiri[2]*rayDiri[2];
	int B = ocVec[0]*rayDiri[0]+ocVec[1]*rayDiri[1]+ocVec[2]*rayDiri[2];
	int C = ocVec[0]*ocVec[0]+ocVec[1]*ocVec[1]+ocVec[2]*ocVec[2];
	A >>= 5; B >>= 4; C >>= 3;

	int disc = B * B - A* (C - ((R*R) >> 3));
	if (disc < 0)  return t;
	int distSqrt = intsqrt(disc);
	int q = (B < 0) ? (-B - distSqrt) : (-B + distSqrt);
	q >>= 1;
	int t0 = (!A) ? 1 << 30 : (q << 12) / A;
	int t1 = (!q) ? 1 << 30 : (C << 10) / q;

	if (t0 > t1)
	{
		int temp = t0; t0 = t1; t1 = temp;
	}
	// if t1 is less than zero, the object is in the ray's negative direction
	// and consequently the ray misses the sphere
	if (t1 < 0) return t1;

	// if t0 is less than zero, the intersection point is at t1
	if (t0 < 0) t = t1;
	else 		t = t0;
	return t;
}

int intabs(int n)
{
	return (n > 0) ? n : -n; 
}

int planeIntersection(int* rayOrg, int* rayDir)
{
	return  -(((65536L*4L-rayOrg[1]) * 4096)/rayDir[1]);
}

// Draw a checker board, position changes over time
u8 checkerBoard(int t, int x, int z)
{
	if (intabs(z) > 65536L*32L) return 0;
	x <<= 4; z <<= 4;
	if (t>0) { x += posX; z += posY; }
	else 	 { x -= posX; z -= posY; }
	x >>= 20; z >>= 20;

	int logoPix = 0;
	if (t > 0) 
	{
		if ((x & 15) >= 4 && (x & 15) < 12)
			logoPix = logo[15-z & 15] & (1 << ((15-x+4) & 7));
	} else
		logoPix = pesthorn[31-(z & 15)*2-(x & 8)/8] & (1 << (x & 7));
	
	if (logoPix) logoPix = 1;
	int checker = ((x >> 4) & 1) ^ ((z >> 4) & 1);
	if (checker)
		logoPix = 1 - logoPix;

	return 128*logoPix;
}

////////////////////////////////////
// Scene 1: Floor and ceiling with a checker board texture
// 	 		and a sphere in between reflecting them
void drawScene1(int lineIdx)
{
	int x,y = lineIdx;
	Veci3 rayOrg, spherePos;
	spherePos[1] = 0;
	spherePos[2] = 2*65536+256*4*intcos(timer/8)+1;
	spherePos[0] = 256*intsin(timer/2);

	int rayDir_y = 2*((y*65536)/YSIZE*XSIZE/YSIZE/2 - 24000);
	int rayOrg_z = -65536L*9/2;
	Veci3 rayDir;

	static int sphereIntersectBuf[XSIZE/4+1];
	for (x = 0; x <= XSIZE/4; x++)
	{
		int rayDir_x = 2*((4*x*65536)/XSIZE - 32768);
		rayDir[0] = (intcos(angle)*rayDir_x - intsin(angle)*rayDir_y ) >> 8;
		rayDir[1] = (intsin(angle)*rayDir_x + intcos(angle)*rayDir_y ) >> 8;
		rayOrg[0] = 0; rayOrg[1] = 0; rayOrg[2] = rayOrg_z;
		rayDir[2] = 65536L;
		sphereIntersectBuf[x] = sphereIntersection(spherePos,rayOrg,rayDir,512);
	}

	int t = 0;
	for (x = 0; x < XSIZE; x++)
	{
		u8 pix = 0;
		rayOrg[0] = 0; rayOrg[1] = 0; rayOrg[2] = rayOrg_z;
		int rayDir_x = 2*((x*65536)/XSIZE - 32768);
		rayDir[0] = (intcos(angle)*rayDir_x - intsin(angle)*rayDir_y ) >> 8;
		rayDir[1] = (intsin(angle)*rayDir_x + intcos(angle)*rayDir_y ) >> 8;
		rayDir[2] = 65536L;

		// If x % 4 == 0, t for sphere intersection was already calculated
		if ((x & 3) == 0) 
			t = sphereIntersectBuf[x/4];
		else
		{
			// No sphere intersection between pixels
			if (sphereIntersectBuf[x/4] < 0 && sphereIntersectBuf[x/4+1] < 0)
				t = -1000;
			// Sphere intersection on both pixels, simple interpolate
			if (sphereIntersectBuf[x/4] > 0 && sphereIntersectBuf[x/4+1] > 0)
				t = sphereIntersectBuf[x/4] + (sphereIntersectBuf[x/4+1]-sphereIntersectBuf[x/4])*(x & 3)/4;
			// Sphere intersection on a single pixel, do an intersection test 
			else
				t = sphereIntersection(spherePos,rayOrg,rayDir,512);
		}

		if (t < 0 || t > 32768) 
		{
			t = (intabs(rayDir[1]) < 8) ? 65536*16 : ((4 << 16) << 12)/rayDir[1];

			if (intabs(t) >= 65536*3/2)
			{
				int plane_x = (rayDir[0]/64 + timer)/32,
				plane_y = (rayDir[1]-8000)/2048+8;
				plane_x = intabs((plane_x+65536) % XSIZE);

				if (plane_y >= 0 && plane_y <= 8)
				  if (nicknameBuf[plane_y*XSIZE/8+plane_x/8] & (1 << (plane_x & 7)))
						pix = 128; else pix = 0;
			} else
			pix = checkerBoard(t,(t >> 8)*(rayDir[0]) >> 4,
					rayOrg_z+((t >> 8)*(rayDir[2]) >> 4)); 
		} else
		{
			rayOrg[0] = rayOrg[0]+ (t*rayDir[0] >> 12);
			rayOrg[1] = rayOrg[1]+ (t*rayDir[1] >> 12);
			rayOrg[2] = rayOrg[2]+ (t*rayDir[2] >> 12);

			Veci3 N;
			N[0] = (rayOrg[0] - spherePos[0]) >> 8;
			N[1] = (rayOrg[1] - spherePos[1]) >> 8;
			N[2] = (rayOrg[2] - spherePos[2]) >> 8;

			int reflProd = (N[0]*rayDir[0] + N[1]*rayDir[1] + N[2]*rayDir[2]) >> 9;
			rayDir[0] = rayDir[0] - (reflProd * N[0] >> 8);
			rayDir[1] = rayDir[1] - (reflProd * N[1] >> 8);
			rayDir[2] = rayDir[2] - (reflProd * N[2] >> 8);

			t = 65536*16;
			if (intabs(rayDir[1]) > 4) 
				t = planeIntersection(rayOrg,rayDir);
			pix = checkerBoard(-t,(t >> 8)*(rayDir[0]) >> 4,
					rayOrg_z+((t >> 8)*(rayDir[2]) >> 4)); 
		}
		int y_byte = (YSIZE-(y+1)) / 8 * XSIZE+(XSIZE-(x+1));
		int y_off = (YSIZE-(y+1)) & 7;

		if (pix > 64)   lcdBuffer[y_byte] |= (1 << y_off);
		else  			lcdBuffer[y_byte] &= ~(1 << y_off);
	}
}


void main_rtr0ket (void) 
{
	int i,x,y;
	for (y = 0; y < 8; y++)
		for (x = 0; x < XSIZE; x++)
			lcdSetPixel(x,y,0);
	DoString(0,0,nickname);

	for (i=0; i<8*XSIZE/8; i++) nicknameBuf[i] = 0; 

	for (y = 0; y < 8; y++)
		for (x = 0; x < XSIZE; x++)
		{
			int x_off = x & 7;
			int pos   = y*XSIZE/8+x/8;
			if (lcdGetPixel(x,y))
				nicknameBuf[pos] |= (1 << x_off);
			else
				nicknameBuf[pos] &= ~(1 << x_off);
		}

	while (1) 
	{
		int key = getInputRaw();
		switch (key) 
		{
			case BTN_ENTER:
				return;
			case BTN_LEFT:
				if (speed < 256) speed += 4;
				break;
			case BTN_RIGHT:
				if (speed > -256) speed -= 4;
				break;
		}
		int dY = 256*256, dX = 256*intsin(timer/32);
		posX += dX*4;
		posY += dY*4;	
		posX &= (1 << 28) - 1;
		posY &= (1 << 28) - 1;
		angle += speed*1000/24/256;
		angle &= 1023;

		int y;
		for (y = 0; y < YSIZE; y++)
			drawScene1(y);	

		lcdDisplay();
		timer += 1000/24;
		delayms_queue_plus(10,0);
	}
}
