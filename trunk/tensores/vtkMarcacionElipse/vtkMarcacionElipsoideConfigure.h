/* 
 * Here is where system computed values get stored.
 * These values should only change when the target compile platform changes.
 */

#define VTKMARCACIONELIPSOIDE_BUILD_SHARED_LIBS
#ifndef VTKMARCACIONELIPSOIDE_BUILD_SHARED_LIBS
#define VTKMARCACIONELIPSOIDE_STATIC
#endif

#if defined(WIN32) && !defined(VTKMARCACIONELIPSOIDE_STATIC)
#pragma warning ( disable : 4275 )

#if defined(vtkMarcacionElipsoide_EXPORTS)
#define VTK_MARCACIONELIPSOIDE_EXPORT __declspec( dllexport ) 
#else
#define VTK_MARCACIONELIPSOIDE_EXPORT __declspec( dllimport ) 
#endif
#else
#define VTK_MARCACIONELIPSOIDE_EXPORT
#endif
