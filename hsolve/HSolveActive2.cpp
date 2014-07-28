/**********************************************************************
** This program is part of 'MOOSE', the
** Messaging Object Oriented Simulation Environment.
**   copyright (C) 2003-2007 Upinder S. Bhalla, Niraj Dudani and NCBS
** It is made available under the terms of the
** GNU Lesser General Public License version 2.1
** See the file COPYING.LIB for the full notice.
**********************************************************************/

#include "header.h"
#include <queue>
#include "HSolveStruct.h"
#include "HSolvePassive2.h"
#include "RateLookup.h"
#include "HSolveActive2.h"
#include "HSolve.h"
#include "../biophysics/CompartmentBase.h"
#include "../biophysics/Compartment.h"
#include "../biophysics/CaConcBase.h"
#include "ZombieCaConc.h"

using namespace moose;
//~ #include "ZombieCompartment.h"
//~ #include "ZombieCaConc.h"

const int HSolveActive::INSTANT_X = 1;
const int HSolveActive::INSTANT_Y = 2;
const int HSolveActive::INSTANT_Z = 4;

HSolveActive::HSolveActive()
{
    caAdvance_ = 1;

    // Default lookup table size
    //~ vDiv_ = 3000;    // for voltage
    //~ caDiv_ = 3000;   // for calcium
}

//////////////////////////////////////////////////////////////////////
// Solving differential equations
//////////////////////////////////////////////////////////////////////
void HSolveActive::step( ProcPtr info )
{
    if ( nCompt_ <= 0 )
        return;

    if ( !current_.size() )
    {
        current_.resize( channel_.size() );
    }

    advanceChannels( info->dt );
    calculateChannelCurrents();
    updateMatrix();
    
    // TODO FastMxElim calculations
    // TODO is build needed?
    FastMatrixElim::advance(V_, passiveOps_, passiveDiagVal_);
    
    advanceCalcium();
    advanceSynChans( info );

    sendValues( info );
    sendSpikes( info );

    externalCurrent_.assign( externalCurrent_.size(), 0.0 );
}

void HSolveActive::calculateChannelCurrents()
{
    vector< ChannelStruct >::iterator ichan;
    vector< CurrentStruct >::iterator icurrent = current_.begin();

    if ( state_.size() != 0 )
    {
        double* istate = &state_[ 0 ];

        for ( ichan = channel_.begin(); ichan != channel_.end(); ++ichan )
        {
            ichan->process( istate, *icurrent );
            ++icurrent;
        }
    }
}

void HSolveActive::updateMatrix() // TODO is it ready?
{
    // resize passiveElim_ // TODO is it needed at every update?
    passiveElim_.setSize(nCompt_, nCompt_); // TODO nCompt_ + 1 cols ???
    
    // copy diagonal and off-diagonal values to be able to alter them
    vector< double > diagvalsCopy = diagvals_;
    map< pair< unsigned int, unsigned int >, double > junctionsCopy = junctions_;
    vector< double > B(nCompt_, .0); // TODO FIXME how is it called ???
    
    // add GkSum and GkEkSum, fill externalCurrent vector
    vector< CurrentStruct >::iterator icurrent = current_.begin();
    vector< currentVecIter >::iterator iboundary = currentBoundary_.begin();
    for (unsigned int i = 0; i < nCompt_; ++i)
    {
        double GkSum = .0, GkEkSum = .0;
        for (; icurrent < *iboundary; ++icurrent)
        {
            GkSum   += icurrent->Gk;
            GkEkSum += icurrent->Gk * icurrent->Ek;
        }
        
        diagvalsCopy[i] += GkSum;
        B[i] = V_[i] * compartment_[i].CmByDt * compartment_[i].EmByRm + GkEkSum;
        
        ++iboundary;
    }
    
    // current injections from inject_
    map< unsigned int, InjectStruct >::iterator inject;
    for (inject = inject_.begin(); inject != inject_.end(); ++inject)
    {
        B[inject->first] +=
            inject->second.injectVarying + inject->second.injectBasal;

        inject->second.injectVarying = 0.0;
    }
    
    // current injections from external sources not included in Hines solver
    for (unsigned int i = 0; i < nCompt_; ++i)
    {
        diagvalsCopy[i] += externalCurrent_[2 * i];
        B[i] += externalCurrent_[2 * i + 1];
    }
    
    // update diagonal values of sparse matrix
    for (unsigned int i = 0; i < nCompt_; ++i)
        passiveElim_.set(i, i, diagvalsCopy[i]);
    // store Gij off-diagonal values into sparse matrix
    map< pair< unsigned int, unsigned int >, double >::iterator jIt;
    for (jIt = junctions_.begin(); jIt != junctions_.end(); ++jIt)
    {
        // TODO is it needed to add both?
        passiveElim_.set(jIt->first.first, jIt->first.second, jIt->second);
        passiveElim_.set(jIt->first.second, jIt->first.first, jIt->second);
    }
    // TODO where to store B ??? the last col of the HS_
}

void HSolveActive::advanceCalcium()
{
    vector< double* >::iterator icatarget = caTarget_.begin();
    vector< double >::iterator ivmid = VMid_.begin();
    vector< CurrentStruct >::iterator icurrent = current_.begin();
    vector< currentVecIter >::iterator iboundary = currentBoundary_.begin();

    caActivation_.assign( caActivation_.size(), 0.0 );

    /*
     * caAdvance_: This flag determines how current flowing into a calcium pool
     * is computed. A value of 0 means that the membrane potential at the
     * beginning of the time-step is used for the calculation. This is how
     * GENESIS does its computations. A value of 1 means the membrane potential
     * at the middle of the time-step is used. This is the correct way of
     * integration, and is the default way.
     */
    if ( caAdvance_ == 1 )
    {
        for ( ; iboundary != currentBoundary_.end(); ++iboundary )
        {
            for ( ; icurrent < *iboundary; ++icurrent )
            {
                if ( *icatarget )
                    **icatarget += icurrent->Gk * ( icurrent->Ek - *ivmid );

                ++icatarget;
            }

            ++ivmid;
        }
    }
    else if ( caAdvance_ == 0 )
    {
        vector< double >::iterator iv = V_.begin();
        double v0;

        for ( ; iboundary != currentBoundary_.end(); ++iboundary )
        {
            for ( ; icurrent < *iboundary; ++icurrent )
            {
                if ( *icatarget )
                {
                    v0 = ( 2 * *ivmid - *iv );

                    **icatarget += icurrent->Gk * ( icurrent->Ek - v0 );
                }

                ++icatarget;
            }

            ++ivmid, ++iv;
        }
    }

    vector< CaConcStruct >::iterator icaconc;
    vector< double >::iterator icaactivation = caActivation_.begin();
    vector< double >::iterator ica = ca_.begin();
    for ( icaconc = caConc_.begin(); icaconc != caConc_.end(); ++icaconc )
    {
        *ica = icaconc->process( *icaactivation );
        ++ica, ++icaactivation;
    }
}

void HSolveActive::advanceChannels( double dt )
{
    vector< double >::iterator iv;
    vector< double >::iterator istate = state_.begin();
    vector< int >::iterator ichannelcount = channelCount_.begin();
    vector< ChannelStruct >::iterator ichan = channel_.begin();
    vector< ChannelStruct >::iterator chanBoundary;
    vector< unsigned int >::iterator icacount = caCount_.begin();
    vector< double >::iterator ica = ca_.begin();
    vector< double >::iterator caBoundary;
    vector< LookupColumn >::iterator icolumn = column_.begin();
    vector< LookupRow >::iterator icarowcompt;
    vector< LookupRow* >::iterator icarow = caRow_.begin();

    LookupRow vRow;
    double C1, C2;
    for ( iv = V_.begin(); iv != V_.end(); ++iv )
    {
        vTable_.row( *iv, vRow );
        icarowcompt = caRowCompt_.begin();
        caBoundary = ica + *icacount;
        for ( ; ica < caBoundary; ++ica )
        {
            caTable_.row( *ica, *icarowcompt );
            ++icarowcompt;
        }

        /*
         * Optimize by moving "if ( instant )" outside the loop, because it is
         * rarely used. May also be able to avoid "if ( power )".
         *
         * Or not: excellent branch predictors these days.
         *
         * Will be nice to test these optimizations.
         */
        chanBoundary = ichan + *ichannelcount;
        for ( ; ichan < chanBoundary; ++ichan )
        {
            if ( ichan->Xpower_ > 0.0 )
            {
                vTable_.lookup( *icolumn, vRow, C1, C2 );
                //~ *istate = *istate * C1 + C2;
                //~ *istate = ( C1 + ( 2 - C2 ) * *istate ) / C2;
                if ( ichan->instant_ & INSTANT_X )
                    *istate = C1 / C2;
                else
                {
                    double temp = 1.0 + dt / 2.0 * C2;
                    *istate = ( *istate * ( 2.0 - temp ) + dt * C1 ) / temp;
                }

                ++icolumn, ++istate;
            }

            if ( ichan->Ypower_ > 0.0 )
            {
                vTable_.lookup( *icolumn, vRow, C1, C2 );
                //~ *istate = *istate * C1 + C2;
                //~ *istate = ( C1 + ( 2 - C2 ) * *istate ) / C2;
                if ( ichan->instant_ & INSTANT_Y )
                    *istate = C1 / C2;
                else
                {
                    double temp = 1.0 + dt / 2.0 * C2;
                    *istate = ( *istate * ( 2.0 - temp ) + dt * C1 ) / temp;
                }

                ++icolumn, ++istate;
            }

            if ( ichan->Zpower_ > 0.0 )
            {
                LookupRow* caRow = *icarow;
                if ( caRow )
                {
                    caTable_.lookup( *icolumn, *caRow, C1, C2 );
                }
                else
                {
                    vTable_.lookup( *icolumn, vRow, C1, C2 );
                }

                //~ *istate = *istate * C1 + C2;
                //~ *istate = ( C1 + ( 2 - C2 ) * *istate ) / C2;
                if ( ichan->instant_ & INSTANT_Z )
                    *istate = C1 / C2;
                else
                {
                    double temp = 1.0 + dt / 2.0 * C2;
                    *istate = ( *istate * ( 2.0 - temp ) + dt * C1 ) / temp;
                }

                ++icolumn, ++istate, ++icarow;
            }
        }

        ++ichannelcount, ++icacount;
    }
}

/**
 * SynChans are currently not under solver's control
 */
void HSolveActive::advanceSynChans( ProcPtr info )
{
    return;
}

void HSolveActive::sendSpikes( ProcPtr info )
{
    vector< SpikeGenStruct >::iterator ispike;
    for ( ispike = spikegen_.begin(); ispike != spikegen_.end(); ++ispike )
        ispike->send( info );
}

/**
 * This function dispatches state values via any source messages on biophysical
 * objects which have been taken over.
 *
 */
void HSolveActive::sendValues( ProcPtr info )
{
    vector< unsigned int >::iterator i;

    for ( i = outVm_.begin(); i != outVm_.end(); ++i )
        moose::Compartment::VmOut()->send(
            //~ ZombieCompartment::VmOut()->send(
            compartmentId_[ *i ].eref(),
            V_[ *i ]
        );

    for ( i = outCa_.begin(); i != outCa_.end(); ++i )
        //~ CaConc::concOut()->send(
        CaConcBase::concOut()->send(
            caConcId_[ *i ].eref(),
            ca_[ *i ]
        );
}
