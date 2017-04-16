#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void ptinpoly1(void *, void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(areapl)(void *, void *, void *, void *);
extern void F77_NAME(inpip)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(k12hat)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kern3d)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(khvc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(khvmat)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(krnqrt)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mse2d)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(n2dist)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nndisf)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nndisg)(void *, void *, void *, void *);
extern void F77_NAME(stkhat)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(stsecal)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(trblik)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(tribble)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(trykh)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"ptinpoly1", (DL_FUNC) &ptinpoly1, 8},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"areapl",  (DL_FUNC) &F77_NAME(areapl),   4},
    {"inpip",   (DL_FUNC) &F77_NAME(inpip),    7},
    {"k12hat",  (DL_FUNC) &F77_NAME(k12hat),  13},
    {"kern3d",  (DL_FUNC) &F77_NAME(kern3d),  13},
    {"khvc",    (DL_FUNC) &F77_NAME(khvc),    16},
    {"khvmat",  (DL_FUNC) &F77_NAME(khvmat),  13},
    {"krnqrt",  (DL_FUNC) &F77_NAME(krnqrt),  16},
    {"mse2d",   (DL_FUNC) &F77_NAME(mse2d),   11},
    {"n2dist",  (DL_FUNC) &F77_NAME(n2dist),   8},
    {"nndisf",  (DL_FUNC) &F77_NAME(nndisf),   7},
    {"nndisg",  (DL_FUNC) &F77_NAME(nndisg),   4},
    {"stkhat",  (DL_FUNC) &F77_NAME(stkhat),  16},
    {"stsecal", (DL_FUNC) &F77_NAME(stsecal), 19},
    {"trblik",  (DL_FUNC) &F77_NAME(trblik),   9},
    {"tribble", (DL_FUNC) &F77_NAME(tribble), 14},
    {"trykh",   (DL_FUNC) &F77_NAME(trykh),   12},
    {NULL, NULL, 0}
};

void R_init_splancs(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
