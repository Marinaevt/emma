//////////////////////////////////////////////////////////////////////
// OpenGLDevice.cpp: Implementierung der Klasse OpenGLDevice.
//
// Copyright by Andrй Stein
// E-Mail: stonemaster@steinsoft.net, andre_stein@web.de
// http://www.steinsoft.net
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
//      Klassenname: OpenGLDevice
//      
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "OpenGLDevice.h"

//////////////////////////////////////////////////////////////////////
// Konstruktion/Destruktion
//////////////////////////////////////////////////////////////////////

OpenGLDevice::OpenGLDevice()
{
	renderContext = nullptr;
	deviceContext = nullptr;
}

OpenGLDevice::~OpenGLDevice()
{
	destroy();
}

OpenGLDevice::OpenGLDevice(HWND& window,int stencil)
{
	create(window,stencil);
}

OpenGLDevice::OpenGLDevice(HDC& deviceContext,int stencil)
{
	create(deviceContext,stencil);
}

//////////////////////////////////////////////////////////////////////
// Цffentliche Funktionen
//////////////////////////////////////////////////////////////////////

bool OpenGLDevice::create(HWND& window,int stencil)
{
	HDC deviceContext_temp = ::GetDC(window);
	
	if (!create(deviceContext_temp,stencil))
	{
		::ReleaseDC(window, deviceContext_temp);

		return false;
	}

	::ReleaseDC(window, deviceContext_temp);
	
	return true;
}

bool OpenGLDevice::create(HDC& deviceContext_temp,int stencil)
{
	if (!deviceContext_temp)
	{
		return false;
	}

	if (!setDCPixelFormat(deviceContext_temp,stencil))
	{
		return false;
	}

	renderContext = wglCreateContext(deviceContext_temp);
	wglMakeCurrent(deviceContext_temp, renderContext);

	OpenGLDevice::deviceContext = deviceContext_temp;
	
	return true;
}



void OpenGLDevice::destroy()
{
	if (renderContext != nullptr)
	{
		wglMakeCurrent(nullptr,nullptr);
		wglDeleteContext(renderContext);
	}
}

void OpenGLDevice::makeCurrent(bool disableOther)
{
	if (renderContext != nullptr)
	{
			//should all other device contexts
			//be disabled then?
			if (disableOther) 
				wglMakeCurrent(nullptr,nullptr);

			wglMakeCurrent(deviceContext, renderContext);
	}
}

//////////////////////////////////////////////////////////////////////
// Protected-Funktionen
//////////////////////////////////////////////////////////////////////


bool OpenGLDevice::setDCPixelFormat(HDC& deviceContext_temp,int stencil)
{
	int pixelFormat;
	DEVMODE resolution;
	//char text[256]; не используется

	//PIXELFORMAT->BPP = DESKTOP->BPP
	EnumDisplaySettings(nullptr, ENUM_CURRENT_SETTINGS, &resolution);
	
	unsigned long dwLayerMask = 0,	// Тоже, что и DWORD в MFC
		dwVisibleMask = 0,
		dwDamageMask = 0;

	static PIXELFORMATDESCRIPTOR pixelFormatDesc =
	{
        sizeof (PIXELFORMATDESCRIPTOR),             
        1,                                          
        PFD_DRAW_TO_WINDOW |                        
        PFD_SUPPORT_OPENGL |
        PFD_DOUBLEBUFFER,
        PFD_TYPE_RGBA,                              
        static_cast<unsigned char>(resolution.dmBitsPerPel),
        0, 0, 0, 0, 0, 0,                           
        0, 
		0,                                       
        0, 
		0, 0, 0, 0,                              
        16,                                         
        static_cast<unsigned char>(stencil),
        0,                                          
        0,                             
        0,                                          
		dwLayerMask, dwVisibleMask, dwDamageMask
    };

    
    
    pixelFormat = ChoosePixelFormat (deviceContext_temp,
									&pixelFormatDesc);

    if (!SetPixelFormat(deviceContext_temp, pixelFormat, 
		&pixelFormatDesc)) 
	{
      return false ;
    }

    return true;
}