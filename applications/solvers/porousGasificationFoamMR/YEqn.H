// eqZx2uHGn041

tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi)")
    )
);

{
    reaction->correct();
    volScalarField Yt(0.0*Y[0]);

    // Species diffusion weighted by porosityF^2 (2026-09-08). With porosityF
    // to the power 1 in both the storage term ddt(porosityF*rho,Yi) and the
    // diffusion term laplacian(porosityF*muEff,Yi), the factor cancels in a
    // region of uniform porosity: the local relaxation time depends on
    // muEff/rho only, whatever the porosity. A nearly solid block (porosityF
    // = 1e-3) was then filled by diffusion within 50-150 s (l^2/D for 2-5 cm
    // and D ~ 3e-5 m2/s).
    //
    // The storage term keeps porosityF^1 (gas inventory in the pores); the
    // diffusion uses porosityF^2, i.e. D_eff = D*eps/tortuosity with a
    // tortuosity ~ 1/eps, a usual engineering approximation of the same
    // family as Bruggeman (eps^1.5) or Millington-Quirk (eps^(4/3)). The
    // original code (porosityF^1 everywhere) amounts to a tortuosity of 1.
    // In the bed (porosityF 0.4-0.69) this slows the internal mixing a little.
    //
    // porosityFf is interpolated with a harmonic mean (resistances in series
    // at the faces).
    surfaceScalarField porosityFf
    (
        "porosityFf",
        1.0/fvc::interpolate(1.0/porosityF)
    );
    surfaceScalarField porosityFf2("porosityFf2", porosityFf*porosityFf);
    surfaceScalarField muEfff("muEfff", fvc::interpolate(turbulence->muEff()));

    // The refractory is a separate region (solidZone) without species
    // equations, so no special treatment of species diffusion is needed at
    // the refractory interface (2026-09-09).

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                (conservativeDdt ? fvm::ddt(porosityF, rho, Yi) : fvm::ddt(porosityF*rho, Yi))
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(porosityFf2 * muEfff, Yi)
             ==
                // 2026-09-28: gas-phase species source weighted by the porosity,
                // like the heat release (EEqn: reaction->Qdot()*porosityF). RR_
                // (laminarEDM) is a rate per unit GAS volume; per unit cell volume
                // it is porosityF*RR. Weighting the heat but not the species broke
                // the energy conservation in the bed (factor eps).
                // Cf. Zuk et al. 2022 (CPC 108407) eq. 2a/2c : meme ponderation.
                porosityF()*reaction->R(Yi)
              + fvOptions(rho, Yi)
              + pyrolysisZone.Srho(i)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            dbgStop("D_YEqn_matrix_" + Yi.name(), dbgMatrix("D_YEqn_matrix_" + Yi.name(), YiEqn));   // NaN check

            YiEqn.solve("Yi");

            dbgStop("D_YEqn_" + Yi.name(), dbgField("D_YEqn_" + Yi.name(), Yi, -0.1, 1.1));   // NaN check

            fvOptions.correct(Yi);

            Yi.max(0.0);
            if (gMax(Yi) > 2.0)
                Info << "warning: mass fraction exceeded 2" << endl;
            Yt += Yi;
        }
        else
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                (conservativeDdt ? fvm::ddt(porosityF, rho, Yi) : fvm::ddt(porosityF*rho, Yi))
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(porosityFf2 * muEfff, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            dbgStop("D_YEqn_matrix_" + Yi.name(), dbgMatrix("D_YEqn_matrix_" + Yi.name(), YiEqn));   // NaN check

            YiEqn.solve("Yi");

            dbgStop("D_YEqn_" + Yi.name(), dbgField("D_YEqn_" + Yi.name(), Yi, -0.1, 1.1));   // NaN check

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    dbgField("D_Yt", Yt, 1e-6, 2.0);   // NaN check

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            Info<< " gas  " << Y[i].name() << " Srho min/max   = "
                << min(pyrolysisZone.Srho(i)()).value() << ", "
                << max(pyrolysisZone.Srho(i)()).value()
                << " min Y = " << min(Y[i]).value() <<" max Y = " << max(Y[i]).value() << endl;

            Y[i] = Y[i] / Yt;

            Info<< " min Ynorm = " << min(Y[i]).value() <<" max Ynorm = " << max(Y[i]).value() << endl;
        }
        else
        {
            Info<< " gas  " <<Y[i].name() << " is inert specie, min Y = "
                << min(Y[i]).value() <<" max Y = " << max(Y[i]).value();

            Y[i] = Y[i]/Yt;

            Info<< " min Ynorm = " << min(Y[i]).value() <<" max Ynorm = " << max(Y[i]).value() << endl;
        }
    }

    // NaN check
    {
        label nYn = 0;
        forAll(Y, i) { nYn += dbgField("D_Ynorm", Y[i], -0.1, 1.1); }
        dbgStop("D_Ynorm", nYn);
    }
}

