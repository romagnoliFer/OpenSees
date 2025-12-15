/* ****************************************************************** **
**   OpenSees - RambergOsgoodSteel (Masing cyclic version, 2025)     **
** ****************************************************************** */

#include <elementAPI.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include "RambergOsgoodSteel.h"

static int numRambergOsgoodSteel = 0;

void* OPS_RambergOsgoodSteel(void)
{
    if (numRambergOsgoodSteel == 0) {
        opserr << "RambergOsgoodSteel (Masing cyclic version) - Stable hysteretic model.
";
        numRambergOsgoodSteel++;
    }

    int iData[1];
    double dData[5]; // Fy, E0, A, N, kMasing
    int numData;

    numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid uniaxialMaterial RambergOsgoodSteel tag
";
        return 0;
    }

    numData = 4;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "WARNING invalid Fy, E0, A, N input
";
        return 0;
    }

    // Optional Masing parameter
    double kMasing = 1.0; // default full Masing
    if (OPS_GetNumRemainingInputArgs() > 0) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData, &kMasing) != 0)
            kMasing = 1.0;
    }

    UniaxialMaterial* theMaterial = new RambergOsgoodSteel(
        iData[0], dData[0], dData[1], dData[2], dData[3], kMasing);
    if (theMaterial == 0)
        opserr << "WARNING could not create RambergOsgoodSteel material
";

    return theMaterial;
}


// Constructors
RambergOsgoodSteel::RambergOsgoodSteel(int tag, double _Fy, double _E0, double _A, double _N, double _kMasing)
    : UniaxialMaterial(tag, MAT_TAG_RambergOsgoodSteel),
      Fy(_Fy), E0(_E0), rezaAA(_A), rezaNN(_N), kMasing(_kMasing)
{
    this->revertToStart();
}

RambergOsgoodSteel::RambergOsgoodSteel(void)
    : UniaxialMaterial(0, MAT_TAG_RambergOsgoodSteel),
      Fy(0.0), E0(0.0), rezaAA(0.002), rezaNN(3.0), kMasing(1.0)
{
    this->revertToStart();
}

RambergOsgoodSteel::~RambergOsgoodSteel(void)
{
}


// Trial strain and hysteretic solver
int RambergOsgoodSteel::setTrialStrain(double trialStrain, double strainRate)
{
    eps = trialStrain;
    double deps = eps - epsP;

    // Detect reversal
    if (deps * (epsP - epsr) < 0.0 && fabs(deps) > 1e-16) {
        epsr = epsP;
        sigr = sigP;
    }

    double sigs0 = (Fy > 0.0) ? Fy : 1.0;

    // Apply Masing scaling: strain halved, stress doubled if kMasing=1
    double scaleEps = 1.0 - 0.5 * kMasing;     // 1→0.5;   0→1.0
    double scaleSig = 1.0 + kMasing;           // 1→2.0;   0→1.0

    double deps_rel = scaleEps * fabs(eps - epsr);

    // Newton–Raphson on modified backbone
    double sig_trial = Fy;
    double sig_old;
    double F, dFdSig;
    const double tol = 1e-8;
    const int maxIter = 50;
    int iter = 0;

    do {
        sig_old = sig_trial;
        F = (sig_trial / E0) + rezaAA * pow(fabs(sig_trial) / sigs0, rezaNN) - deps_rel;
        dFdSig = (1.0 / E0) + rezaAA * (rezaNN / sigs0) * pow(fabs(sig_trial) / sigs0, rezaNN - 1.0);
        sig_trial -= F / dFdSig;
        iter++;
    } while (fabs(sig_trial - sig_old) > tol && iter < maxIter);

    if (iter == maxIter)
        opserr << "RambergOsgoodSteel::setTrialStrain -- NR did not converge at strain " << eps << endln;

    double sign = (eps - epsr >= 0.0) ? 1.0 : -1.0;
    sig = sigr + sign * scaleSig * sig_trial;

    e = scaleSig / ((1.0 / E0) + rezaAA * (rezaNN / sigs0) * pow(fabs(sig_trial) / sigs0, rezaNN - 1.0));
    return 0;
}


// Accessors
double RambergOsgoodSteel::getStrain(void) { return eps; }
double RambergOsgoodSteel::getStress(void) { return sig; }
double RambergOsgoodSteel::getTangent(void) { return e; }
double RambergOsgoodSteel::getInitialTangent(void) { return E0; }


// State handling
int RambergOsgoodSteel::commitState(void)
{
    epsP = eps;
    sigP = sig;
    eP = e;
    return 0;
}

int RambergOsgoodSteel::revertToLastCommit(void)
{
    eps = epsP;
    sig = sigP;
    e = eP;
    return 0;
}

int RambergOsgoodSteel::revertToStart(void)
{
    eps = 0.0;
    sig = 0.0;
    e = E0;
    epsP = sigP = 0.0;
    eP = E0;
    epsr = sigr = 0.0;
    return 0;
}


// Copy and print
UniaxialMaterial* RambergOsgoodSteel::getCopy(void)
{
    RambergOsgoodSteel* theCopy = new RambergOsgoodSteel(this->getTag(), Fy, E0, rezaAA, rezaNN, kMasing);
    *theCopy = *this;
    return theCopy;
}

void RambergOsgoodSteel::Print(OPS_Stream& s, int flag)
{
    s << "RambergOsgoodSteel, tag: " << this->getTag() << endln;
    s << "  Fy = " << Fy << ", E0 = " << E0 << ", A = " << rezaAA << ", N = " << rezaNN
      << ", kMasing = " << kMasing << endln;
    s << "  strain: " << eps << ", stress: " << sig << ", tangent: " << e << endln;
}


// Parallel channel support
int RambergOsgoodSteel::sendSelf(int commitTag, Channel& theChannel)
{
    static Vector data(10);
    data(0) = Fy;
    data(1) = E0;
    data(2) = rezaAA;
    data(3) = rezaNN;
    data(4) = kMasing;
    data(5) = epsP;
    data(6) = sigP;
    data(7) = eP;
    data(8) = epsr;
    data(9) = sigr;

    if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "RambergOsgoodSteel::sendSelf -- failed
";
        return -1;
    }
    return 0;
}

int RambergOsgoodSteel::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    static Vector data(10);
    if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
        opserr << "RambergOsgoodSteel::recvSelf -- failed
";
        return -1;
    }

    Fy = data(0);
    E0 = data(1);
    rezaAA = data(2);
    rezaNN = data(3);
    kMasing = data(4);
    epsP = data(5);
    sigP = data(6);
    eP = data(7);
    epsr = data(8);
    sigr = data(9);

    eps = epsP;
    sig = sigP;
    e = eP;
    return 0;
}