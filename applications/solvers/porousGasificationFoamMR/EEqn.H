// eqZx2uHGn044

if (solveEnergy)
{
    if (pyrolysisZone.equilibrium())
    {}
    else
    {
        volScalarField& he = thermo.he();
        volScalarField heatUpGas = pyrolysisZone.heatUpGas()();
        heatTransfer = pyrolysisZone.heatTransfer()();

        // 2026-10-12 (energy balance, paper Section 3.7): sensible enthalpy of the injected mass.
        // With the mass source Srho of rhoEqn, ddt(eps*rho*he) + div(phi*he) adds the mass at he = 0
        // (298 K). Produced species (Srho(i) > 0) enter at Ts, consumed species (O2, Srho(i) < 0)
        // leave at the gas temperature T.
        volScalarField injSh
        (
            IOobject("injSh", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            mesh,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        );
        {
            const volScalarField& TsSol = pyrolysisZone.T();
            basicSpecieMixture& compoInj = thermo.composition();
            forAll(Y, i)
            {
                const volScalarField SrhoI(pyrolysisZone.Srho(i));
                forAll(injSh, celli)
                {
                    const scalar s = SrhoI[celli];
                    if (s != 0)
                    {
                        const scalar Ti = (s > 0) ? max(TsSol[celli], 200.0) : T[celli];
                        injSh[celli] += s*compoInj.Hs(i, p[celli], Ti);
                    }
                }
            }
        }
        // 2026-10-12: gas emission, not removed by radiation->Sh() (eCont = 0).
        const volScalarField gasEmission
        (
            "gasEmission",
            4.0*radiation->heterogeneousAbsorptionEmission().aCont()*constant::physicoChemical::sigma*pow4(T)
        );
        Info<< " _cons injSh int = " << gSum((injSh*mesh.V())()) << " W ; gasEmission int = "
            << gSum((gasEmission*mesh.V())()) << " W" << endl;

        Info<< " rho gas min/max = " << gMin(rho) << ", " << gMax(rho) << endl;
        Info<< " Sh gas min/max = " << gMin(chemistrySh) << ", " << gMax(chemistrySh) << endl;
        Info<< " porosity min/max = " << gMin(porosityF) << ", " << gMax(porosityF) << endl;
        Info<< " heatUpGas min/max = " << gMin(heatUpGas) << ", " << gMax(heatUpGas) << endl;
        Info<< " heatTransfer min/max = " << gMin(heatTransfer) << ", " << gMax(heatTransfer) << endl;

        // porosityF on divq(he) (2026-09-03): the thermal diffusion was not
        // weighted by the porosity, unlike the ddt and convection terms.
        // (The same weighting on heatTransfer was removed on 2026-10-12, see
        // below: it is one of the defects described in the paper, Section 3.7.)
        fvScalarMatrix EEqn
        (
            fvm::ddt(porosityF * rho, he)
          + fvm::div(phi, he)
          + fvc::ddt(porosityF * rho, K)
          + fvc::div(phi, K)
          + (
                he.name() == "e"
              ? fvc::div
                (
                    fvc::absolute(phi/fvc::interpolate(porosityF * rho), U),
                    p,
                    "div(phiv,p)"
                )
              : -dpdt
            )
          + porosityF * thermophysicalTransport->divq(he)
         ==
            reaction->Qdot() * porosityF
          + fvOptions(rho, he)
          + heatTransfer          // 2026-10-12: without porosityF (upstream, conservative form)
          + heatUpGas             // zero (heatUpGasCalc returns 0 in this library), replaced by injSh
          + injSh                 // pyrolysis gas injected with its sensible enthalpy
          + radiation->Sh(thermo, he)
          - gasEmission           // gas emission 4 a sigma T^4 removed (eCont = 0 in the model)
        );

        EEqn.relax();

        fvOptions.constrain(EEqn);

        // NaN check, stage E: source terms and matrix before the solution
        {
            label nE = 0;
            {
                const volScalarField QdotDbgE("QdotE", reaction->Qdot());
                nE += dbgField("E_src", QdotDbgE, -1e10, 1e10);
            }
            nE += dbgField("E_src", heatTransfer, -1e9, 1e9);
            nE += dbgField("E_src", injSh, -1e10, 1e10);
            nE += dbgField("E_src", gasEmission, -1e10, 1e10);
            nE += dbgField("E_src", he, -1e7, 1e8);
            dbgStop("E_src", nE);
            dbgStop("E_matrix", dbgMatrix("E_matrix", EEqn));
        }

        EEqn.solve();

        dbgStop("E_he", dbgField("E_he", he, -1e6, 1e8));   // NaN check

        fvOptions.correct(he);

        thermo.correct();

        dbgStop("E_thermoT", dbgField("E_thermoT", T, 250, 3500));   // NaN check

        // SIGFPE diagnostic of hePsiThermo::calculate (2026-10-05), debug build
        // only; run with FOAM_SIGFPE=false.
        {
            const volScalarField muDbg(thermo.mu());
            const volScalarField alphaDbg(thermo.alpha());
            const volScalarField CpDbg(thermo.Cp());
            const volScalarField QdotDbg(reaction->Qdot());
            label nBad = 0;
            forAll(T, celli)
            {
                const scalar Tc = T[celli];
                const bool bad =
                    !std::isfinite(Tc) || !std::isfinite(he[celli])
                 || !std::isfinite(muDbg[celli]) || !std::isfinite(alphaDbg[celli])
                 || !std::isfinite(CpDbg[celli]) || CpDbg[celli] <= 0
                 || Tc < 250 || Tc > 3500;
                if (bad && nBad < 20)
                {
                    scalar sumY = 0;
                    forAll(Y, i) { sumY += Y[i][celli]; }
                    Pout<< "DBG BAD cell " << celli << " C " << mesh.C()[celli]
                        << " V " << mesh.V()[celli]
                        << " he " << he[celli] << " heOld " << he.oldTime()[celli]
                        << " T " << Tc << " Told " << T.oldTime()[celli]
                        << " Cp " << CpDbg[celli] << " mu " << muDbg[celli]
                        << " alpha " << alphaDbg[celli]
                        << " p " << p[celli] << " rho " << rho[celli]
                        << " eps " << porosityF[celli]
                        << " heatUpGas " << heatUpGas[celli]
                        << " heatTransfer " << heatTransfer[celli]
                        << " Qdot " << QdotDbg[celli]
                        << " sumY " << sumY << nl << "DBG   Y:";
                    forAll(Y, i) { Pout<< " " << Y[i].name() << "=" << Y[i][celli]; }
                    Pout<< nl << "DBG   Yold:";
                    forAll(Y, i) { Pout<< " " << Y[i].name() << "=" << Y[i].oldTime()[celli]; }
                    Pout<< endl;
                }
                if (bad) nBad++;
            }
            reduce(nBad, sumOp<label>());
            if (nBad > 0)
            {
                Info<< "DBG nBad = " << nBad << " -> writeNow + stop" << endl;
                runTime.writeNow();
                FatalErrorInFunction << "DBG stop after bad thermo cells" << exit(FatalError);
            }
        }

        Info<< " hs min/max = " << gMin(he) << ", " << gMax(he) << endl;
        Info<< " T gas min/max = " << gMin(T) << ", "
            << gMax(T) << endl;
    }
}
