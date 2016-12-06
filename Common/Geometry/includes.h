// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once


// Including SDKDDKVer.h defines the highest available Windows platform.

// If you wish to build your application for a previous Windows platform, include WinSDKVer.h and
// set the _WIN32_WINNT macro to the platform you wish to support before including SDKDDKVer.h.

//#include <SDKDDKVer.h>

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>

// TODO: reference additional headers your program requires here

#include <memory>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
#include <tuple>
#include <cassert>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <random>
#include <chrono>