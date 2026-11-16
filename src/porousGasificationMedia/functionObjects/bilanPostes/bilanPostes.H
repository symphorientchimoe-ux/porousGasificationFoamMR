/*---------------------------------------------------------------------------*\
    bilanPostes (16/11/2026) : bilan d'energie par poste, article ECM.

    Integre a chaque pas de temps (somme flux*dt) des flux de frontiere et ecrit
    leur moyenne sur chaque intervalle « interval » (s simulees) dans
    postProcessing/bilanPostes/<temps de depart>/<mode>.dat (processeur maitre).

    mode fluide : t duree phi_in Ha_in phi_out Ha_out Hs_out m_<espece>...
        phi en kg/s (phi_in < 0), H en W (ha : enthalpie absolue, hs : sensible),
        m_<espece> = flux massique sortant de l'espece (kg/s), patchs inlet et outlet.
    mode solide : t duree q_paroi
        q_paroi (W) = somme kappa snGrad(T) |Sf| sur le patch outerWall (< 0 : perte).

    Exemple :
        bilanFluide
        {
            type     bilanPostes;
            libs     ("libbilanPostes.so");
            region   fluidZone;
            mode     fluide;
            species  (CO H2 CH4 C6H5OH CO2 H2O O2);
            interval 1;
        }

    Remplace les functionObjects « coded » du meme nom : OpenFOAM 8 refuse le code
    compile a la volee quand le calcul tourne en root (conteneur des runners GitHub).
\*---------------------------------------------------------------------------*/

#ifndef functionObjects_bilanPostes_H
#define functionObjects_bilanPostes_H

#include "fvMeshFunctionObject.H"
#include "OFstream.H"
#include "scalarList.H"
#include "wordList.H"

namespace Foam
{
namespace functionObjects
{

class bilanPostes
:
    public fvMeshFunctionObject
{
    // Donnees

        word mode_;
        wordList species_;
        word inletPatch_;
        word outletPatch_;
        word wallPatch_;
        scalar interval_;
        scalar duree_;
        scalarList acc_;
        autoPtr<OFstream> os_;

    // Fonctions privees

        void ouvrir();
        void flux(scalarList& f) const;

public:

    TypeName("bilanPostes");

    bilanPostes
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~bilanPostes();

    virtual bool read(const dictionary&);
    virtual bool execute();
    virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif
