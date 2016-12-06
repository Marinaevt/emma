// OpenGLDevice.h: Schnittstelle fьr die Klasse OpenGLDevice.
//
// Copyright by Andrй Stein
// E-Mail: stonemaster@steinsoft.net, andre_stein@web.de
// http://www.steinsoft.net
//////////////////////////////////////////////////////////////////////

#pragma once
#include <windows.h>

class OpenGLDevice
{
public:
	OpenGLDevice(HDC& deviceContext,int stencil = 0);
	OpenGLDevice(HWND& window,int stencil = 0);
	OpenGLDevice();

	bool create(HDC& deviceContext,int  stencil = 0);
	bool create(HWND& window,int stencil = 0);

	void destroy();
	void makeCurrent(bool disableOther = true);

	
	virtual ~OpenGLDevice();

protected:
	bool setDCPixelFormat(HDC& deviceContext,int stencil);
	
	HGLRC renderContext;
	HDC deviceContext;
};

