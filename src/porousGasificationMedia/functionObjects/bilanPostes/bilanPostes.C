/*---------------------------------------------------------------------------*\
    bilanPostes (16/11/2026) : voir bilanPostes.H.
\*---------------------------------------------------------------------------*/

#include "bilanPostes.H"
#include "basicThermo.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(bilanPostes, 0);
    addToRunTimeSelectionTable(functionObject, bilanPostes, dictionary);
}
}


void Foam::functionObjects::bilanPostes::ouvrir()
{
    if (!Pstream::master() || os_.valid())
    {
        return;
    }
    const Time& rt = time_;
    const fileName d
    (
        rt.rootPath()/rt.globalCaseName()/"postProcessing"/"bilanPostes"
       /rt.timeName()
    );
    mkDir(d);
    os_.reset(new OFstream(d/(mode_ + ".dat")));
    os_().precision(12);
    if (mode_ == "fluide")
    {
        os_() << "# t duree phi_in Ha_in phi_out Ha_out Hs_out";
        forAll(species_, i)
        {
            os_() << " m" << species_[i];
        }
        os_() << endl;
    }
    else
    {
        os_() << "# t duree q_paroi" << endl;
    }
}


void Foam::functionObjects::bilanPostes::flux(scalarList& f) const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const basicThermo& th =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    if (mode_ == "fluide")
    {
        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>("phi");
        const label pin = mesh_.boundaryMesh().findPatchID(inletPatch_);
        const label pout = mesh_.boundaryMesh().findPatchID(outletPatch_);
        const scalarField& phin = phi.boundaryField()[pin];
        const scalarField& phout = phi.boundaryField()[pout];
        const scalarField hain(th.ha(T.boundaryField()[pin], pin));
        const scalarField haout(th.ha(T.boundaryField()[pout], pout));
        const scalarField hsout(th.hs(T.boundaryField()[pout], pout));

        f.setSize(5 + species_.size());
        f[0] = gSum(phin);
        f[1] = gSum(phin*hain);
        f[2] = gSum(phout);
        f[3] = gSum(phout*haout);
        f[4] = gSum(phout*hsout);
        forAll(species_, i)
        {
            const volScalarField& Y =
                mesh_.lookupObject<volScalarField>(species_[i]);
            f[5 + i] = gSum(phout*Y.boundaryField()[pout]);
        }
    }
    else
    {
        const label pw = mesh_.boundaryMesh().findPatchID(wallPatch_);
        const tmp<volScalarField> tk(th.kappa());
        f.setSize(1);
        f[0] = gSum
        (
            tk().boundaryField()[pw]*T.boundaryField()[pw].snGrad()
           *mesh_.magSf().boundaryField()[pw]
        );
    }
}


Foam::functionObjects::bilanPostes::bilanPostes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    mode_("fluide"),
    species_(),
    inletPatch_("inlet"),
    outletPatch_("outlet"),
    wallPatch_("outerWall"),
    interval_(1),
    duree_(0),
    acc_()
{
    read(dict);
}


Foam::functionObjects::bilanPostes::~bilanPostes()
{}


bool Foam::functionObjects::bilanPostes::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    mode_ = dict.lookupOrDefault<word>("mode", "fluide");
    species_ = dict.lookupOrDefault<wordList>("species", wordList());
    inletPatch_ = dict.lookupOrDefault<word>("inletPatch", "inlet");
    outletPatch_ = dict.lookupOrDefault<word>("outletPatch", "outlet");
    wallPatch_ = dict.lookupOrDefault<word>("wallPatch", "outerWall");
    interval_ = dict.lookupOrDefault<scalar>("interval", 1.0);
    if (mode_ != "fluide" && mode_ != "solide")
    {
        FatalIOErrorInFunction(dict)
            << "mode " << mode_ << " inconnu (fluide ou solide)"
            << exit(FatalIOError);
    }
    return true;
}


bool Foam::functionObjects::bilanPostes::execute()
{
    ouvrir();

    const scalar dt = time_.deltaTValue();
    scalarList f;
    flux(f);
    if (acc_.size() != f.size())
    {
        acc_.setSize(f.size(), 0.0);
        acc_ = 0.0;
    }
    forAll(f, i)
    {
        acc_[i] += f[i]*dt;
    }
    duree_ += dt;

    if (duree_ >= interval_ - 1e-9)
    {
        if (Pstream::master())
        {
            os_() << time_.value() << " " << duree_;
            forAll(acc_, i)
            {
                os_() << " " << acc_[i]/duree_;
            }
            os_() << endl;
        }
        acc_ = 0.0;
        duree_ = 0;
    }
    return true;
}


bool Foam::functionObjects::bilanPostes::write()
{
    return true;
}
