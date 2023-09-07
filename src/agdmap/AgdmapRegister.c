/**CFile****************************************************************

  FileName    [AgdmapRegister.c]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Agdmap.]

  Synopsis    [External declarations.]

  Author      [Longfei Fan]

  Affiliation [Fudan University]

  Date        [August 1, 2023]

  Vision      [ 1.0 ]

***********************************************************************/

#include "base/main/main.h"
#include "AgdmapCommand.h"

ABC_NAMESPACE_HEADER_START

////////////////////////////////////////////////////////////////////////
///                      DECLARATIONS                                ///
////////////////////////////////////////////////////////////////////////

static int Agdmap_Mapping(Abc_Frame_t * pAbc, int argc, char **argv);

////////////////////////////////////////////////////////////////////////
///                      FUNCTION DEFINITIONS                        ///
////////////////////////////////////////////////////////////////////////

/**Function*************************************************************

  Synopsis    [initialisation procedure]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Agdmap_Init(Abc_Frame_t * pAbc)
{
	Cmd_CommandAdd(pAbc, "FPGA mapping", "agdmap",    Agdmap_Mapping,        1);
}
/**Function*************************************************************

  Synopsis    [initialisation procedure

  Description []
        
  SideEffects []

  SeeAlso     []

***********************************************************************/
static int  Agdmap_Mapping(Abc_Frame_t * pAbc, int argc, char **argv)
{
	return Agdmap(pAbc, argc, argv);
}

////////////////////////////////////////////////////////////////////////
///                              END                                 ///
////////////////////////////////////////////////////////////////////////

ABC_NAMESPACE_HEADER_END
