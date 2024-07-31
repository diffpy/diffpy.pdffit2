/***********************************************************************
* Python 2.7 needs MSVC 9.0 which does not have unique_ptr.
* Work around by substituting auto_ptr instead.
***********************************************************************/

#ifdef _MSC_VER

// workarounds for MSVC 9.0 --------------------------------------------
#if _MSC_VER <= 1500
#define unique_ptr auto_ptr
#endif

#endif  // _MSC_VER
