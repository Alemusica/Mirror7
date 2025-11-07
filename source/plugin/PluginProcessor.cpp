#include "plugin/PluginProcessor.h"
#include "gui/PluginEditor.h"
#include "engine/Mirror7Engine.h"
#include "config/ParamIDs.h"

#include <cmath>

using namespace ParamIDs;

static void mirror7Log (const juce::String& msg)
{
    // Best-effort file log in user AppData
    auto f = juce::File::getSpecialLocation (juce::File::userApplicationDataDirectory)
               .getChildFile ("Mirror7.log");
    if (auto stream = f.createOutputStream())
    {
        const auto line = juce::Time::getCurrentTime().toString (true, true) + ": " + msg + "\n";
        stream->writeText (line, false, false, "UTF-8");
        stream->flush();
    }
}

//==============================================================================
Mirror7AudioProcessor::Mirror7AudioProcessor()
    : parameters (*this, nullptr, juce::Identifier("mirror7"), createParameterLayout())
{
    mirror7Log ("ctor");
    engine = std::make_unique<Mirror7Engine>();
}

void Mirror7AudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    mirror7Log ("prepareToPlay sr=" + juce::String (sampleRate) + " n=" + juce::String (samplesPerBlock));
    if (!engine)
        engine = std::make_unique<Mirror7Engine>();

    const double sr = (sampleRate > 0.0 ? sampleRate : 44100.0);
    const int blockSize = juce::jmax (1, samplesPerBlock);

    engine->prepare (sr, blockSize);
    engineSampleRate = sr;
    engineBlockSize = blockSize;
    enginePrepared = true;
    const int osSel = (int) parameters.getRawParameterValue (osFactor)->load();
    engine->setOversampling (osSel);
}

void Mirror7AudioProcessor::releaseResources()
{
    if (engine)
        engine->reset();

    enginePrepared = false;
    engineSampleRate = 0.0;
    engineBlockSize = 0;
}

bool Mirror7AudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;
    return true;
}

void Mirror7AudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer&)
{
    juce::ScopedNoDenormals noDenormals;
    if (!engine)
        engine = std::make_unique<Mirror7Engine>();

    const double hostSampleRate = getSampleRate();
    const int hostBlockSize = getBlockSize();

    const int currentBlockSamples = buffer.getNumSamples();

    if (!enginePrepared
        || (hostSampleRate > 0.0 && std::abs (hostSampleRate - engineSampleRate) > 1.0e-6)
        || (hostBlockSize > 0 && hostBlockSize != engineBlockSize)
        || (hostBlockSize <= 0 && currentBlockSamples > engineBlockSize))
    {
        const double sr = (hostSampleRate > 0.0 ? hostSampleRate
                                                : (engineSampleRate > 0.0 ? engineSampleRate : 44100.0));
        const int blockSize = juce::jmax (1, hostBlockSize > 0 ? hostBlockSize : currentBlockSamples);
        engine->prepare (sr, blockSize);
        engineSampleRate = sr;
        engineBlockSize = blockSize;
        enginePrepared = true;
    }

    if (!enginePrepared)
    {
        buffer.clear();
        return;
    }

    const float outGainDb = parameters.getRawParameterValue (outGain)->load();
    const float panv      = parameters.getRawParameterValue (pan)->load();
    const float widthv    = parameters.getRawParameterValue (width)->load();
    const float itdUsv    = parameters.getRawParameterValue (itdUs)->load();
    const float ildDbv    = parameters.getRawParameterValue (ildDb)->load();
    const float velAmtv   = parameters.getRawParameterValue (velvetAmt)->load();
    const bool  autoPanOn = parameters.getRawParameterValue (autopan)->load() > 0.5f;
    const float autoSpv   = parameters.getRawParameterValue (autopanSp)->load();
    const int   osSel     = (int) parameters.getRawParameterValue (osFactor)->load();

    if (engine && engine->getOversamplingFactor() != (osSel == 0 ? 1 : (osSel == 1 ? 2 : 4)))
        engine->setOversampling (osSel);

    // Map parameters into Engine API
    engine->setVelvetAmount (velAmtv);
    engine->setPanWidth (panv, widthv);
    engine->setITDILD (itdUsv, ildDbv);
    engine->setAutopan (autoPanOn, autoSpv);

    // Update full Mirror7 params
    auto& P = engine->getParams();
    P.rateHz = parameters.getRawParameterValue (rateHz)->load();
    P.baseLenMs = parameters.getRawParameterValue (baseMs)->load();
    P.lenPhi = parameters.getRawParameterValue (lenPhi)->load();
    P.hemisCoupling = parameters.getRawParameterValue (hemisCoupling)->load();
    P.spatMinDeg = parameters.getRawParameterValue (spatMinDeg)->load();
    P.spatMinMs  = parameters.getRawParameterValue (spatMinMs)->load();
    P.spatIPD    = parameters.getRawParameterValue (spatIPD)->load();
    P.spatShadow = parameters.getRawParameterValue (spatShadow)->load();
    P.spatMirror = parameters.getRawParameterValue (spatMirror)->load() > 0.5f;
    P.envAttack  = parameters.getRawParameterValue (envA)->load();
    P.envDecay   = parameters.getRawParameterValue (envD)->load();
    P.envSustain = parameters.getRawParameterValue (envS)->load();
    P.envRelease = parameters.getRawParameterValue (envR)->load();
    P.dialogueOn = parameters.getRawParameterValue (dialogOn)->load() > 0.5f;
    P.dialogueStrength = parameters.getRawParameterValue (dialogStr)->load();
    P.dialogueMemory   = parameters.getRawParameterValue (dialogMem)->load();
    P.dialoguePhiMix   = parameters.getRawParameterValue (dialogPhi)->load();
    P.glitchMix = parameters.getRawParameterValue (glitchMix)->load();
    P.srCrush   = parameters.getRawParameterValue (srCrush)->load();
    P.bitCrush  = parameters.getRawParameterValue (bitCrush)->load();
    P.phiPan        = parameters.getRawParameterValue (phiPan)->load() > 0.5f;
    P.phiMode       = parameters.getRawParameterValue (phiMode)->load() > 0.5f;
    P.noiseMode     = (int) parameters.getRawParameterValue (noiseMode)->load();
    P.noiseColor    = (int) parameters.getRawParameterValue (noiseColor)->load();
    P.colorAmt      = parameters.getRawParameterValue (colorAmt)->load();

    // Aureo / Quantum noise families
    P.aureoMix      = parameters.getRawParameterValue (aureoMix)->load();
    P.aureoDecay    = parameters.getRawParameterValue (aureoDecay)->load();
    P.aureoStride   = parameters.getRawParameterValue (aureoStride)->load();
    P.aureoHarmonics= parameters.getRawParameterValue (aureoHarmo)->load();
    P.aureoVelvet   = parameters.getRawParameterValue (aureoVelvet)->load() > 0.5f;

    P.qMix          = parameters.getRawParameterValue (qMix)->load();
    P.qDetail       = parameters.getRawParameterValue (qDetail)->load();
    P.qStride       = parameters.getRawParameterValue (qStride)->load();
    P.qBase         = parameters.getRawParameterValue (qBase)->load();
    P.qVelvet       = parameters.getRawParameterValue (qVelvet)->load();
    P.vhsWow        = parameters.getRawParameterValue (vhsWow)->load();
    P.vhsFlutter    = parameters.getRawParameterValue (vhsFlutter)->load();

    // Modal layer
    P.modalOn       = parameters.getRawParameterValue (modalOn)->load() > 0.5f;
    P.modalMirror   = parameters.getRawParameterValue (modalMirror)->load() > 0.5f;
    P.modalMix      = parameters.getRawParameterValue (modalMix)->load();
    P.modalDecay    = parameters.getRawParameterValue (modalDecay)->load();
    P.modalPreset   = (int) parameters.getRawParameterValue (modalPreset)->load();

    // Φ geometry fine controls
    P.phiRatio      = parameters.getRawParameterValue (phiRatio)->load();
    P.phiHeadB      = parameters.getRawParameterValue (phiHeadB)->load();
    P.phiDistance   = parameters.getRawParameterValue (phiDistance)->load();
    P.phiElev       = parameters.getRawParameterValue (phiElev)->load();
    P.phiElevNotch  = parameters.getRawParameterValue (phiElevNotch)->load();
    P.phiTorsoMix   = parameters.getRawParameterValue (phiTorsoMix)->load();
    P.phiTorsoMs    = parameters.getRawParameterValue (phiTorsoMs)->load();
    P.phiTorsoHpHz  = parameters.getRawParameterValue (phiTorsoHp)->load();

    // Spatial profile / pinna / CEFG
    P.pinnaOn       = parameters.getRawParameterValue (pinnaOn)->load() > 0.5f;
    P.pinnaDepth    = parameters.getRawParameterValue (pinnaDepth)->load();
    P.spatialProfile= (int) parameters.getRawParameterValue (spatialProfile)->load();
    P.cefgGain      = parameters.getRawParameterValue (cefgGain)->load();
    P.cefgMix       = parameters.getRawParameterValue (cefgMix)->load();

    // Host sync
    P.syncEnable    = parameters.getRawParameterValue (syncEnable)->load() > 0.5f;
    P.syncDivision  = parameters.getRawParameterValue (syncDivision)->load();
    P.syncSlew      = parameters.getRawParameterValue (syncSlew)->load();

    // Host transport info
    double hostBpm = -1.0; bool isPlaying = false;
    if (auto* ph = getPlayHead())
    {
        if (auto pos = ph->getPosition())
        {
            if (auto bpm = pos->getBpm()) hostBpm = *bpm;
            isPlaying = pos->getIsPlaying();
        }
    }
    engine->updateHostTransport (hostBpm, isPlaying);

    engine->process (buffer, outGainDb);
}

//==============================================================================
juce::AudioProcessorEditor* Mirror7AudioProcessor::createEditor()
{
    return new Mirror7AudioProcessorEditor (*this);
}

//==============================================================================
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new Mirror7AudioProcessor();
}

//==============================================================================
void Mirror7AudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    juce::MemoryOutputStream os (destData, false);
    parameters.state.writeToStream (os);
}

void Mirror7AudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    auto tree = juce::ValueTree::readFromData (data, (size_t) sizeInBytes);
    if (tree.isValid())
        parameters.replaceState (tree);
}

//==============================================================================
juce::AudioProcessorValueTreeState::ParameterLayout Mirror7AudioProcessor::createParameterLayout()
{
    std::vector<std::unique_ptr<juce::RangedAudioParameter>> p;

    p.push_back (std::make_unique<juce::AudioParameterFloat>(outGain, "Output", juce::NormalisableRange<float> (-24.0f, 12.0f, 0.01f), -6.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(pan,     "Pan",    juce::NormalisableRange<float> (-1.0f, 1.0f, 0.0001f), 0.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(width,   "Width",  juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 1.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(itdUs,   "ITD (us)", juce::NormalisableRange<float> (0.0f, 1600.0f, 0.01f), 700.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(ildDb,   "ILD (dB)", juce::NormalisableRange<float> (0.0f, 30.0f, 0.01f), 20.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(velvetAmt,"Velvet", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.12f));
    p.push_back (std::make_unique<juce::AudioParameterBool>(autopan,  "Autopan", true));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(autopanSp,"Autopan Speed", juce::NormalisableRange<float> (0.02f, 2.0f, 0.0001f, 0.3f), 0.12f));
    p.push_back (std::make_unique<juce::AudioParameterChoice>(osFactor, "Oversample", juce::StringArray { "1x", "2x", "4x" }, 0));

    // Mirror7 core params
    p.push_back (std::make_unique<juce::AudioParameterFloat>(rateHz,   "Rate (Hz)", juce::NormalisableRange<float> (0.1f, 240.0f, 0.0001f, 0.3f), 8.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(baseMs,   "Base Len (ms)", juce::NormalisableRange<float> (1.0f, 2000.0f, 0.01f, 0.35f), 120.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(lenPhi,   "Len Φ", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.8f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(hemisCoupling, "Hemis Coupling", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.6f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(spatMinDeg, "Min Deg", juce::NormalisableRange<float> (0.0f, 90.0f, 0.01f), 12.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(spatMinMs,  "Min ms",  juce::NormalisableRange<float> (0.0f, 500.0f, 0.01f), 35.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(spatIPD,    "IPD",      juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.6f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(spatShadow, "Shadow",   juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.7f));
    p.push_back (std::make_unique<juce::AudioParameterBool>(spatMirror,  "Mirror",   false));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(envA,      "Env A",    juce::NormalisableRange<float> (0.01f, 4.0f, 0.0001f, 0.35f), 0.18f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(envD,      "Env D",    juce::NormalisableRange<float> (0.01f, 4.0f, 0.0001f, 0.35f), 0.28f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(envS,      "Env S",    juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.55f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(envR,      "Env R",    juce::NormalisableRange<float> (0.01f, 4.0f, 0.0001f, 0.35f), 0.30f));
    p.push_back (std::make_unique<juce::AudioParameterBool>(dialogOn,   "Dialogue", true));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(dialogStr, "Dial Str", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.6f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(dialogMem, "Dial Mem", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.5f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(dialogPhi, "Dial ΦMix",juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.75f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(glitchMix, "Glitch",   juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.5f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(srCrush,   "SR Crush", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.2f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(bitCrush,  "BitCrush", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.15f));
    p.push_back (std::make_unique<juce::AudioParameterBool>(phiPan,     "Φ Pan",    false));
    p.push_back (std::make_unique<juce::AudioParameterBool>(phiMode,    "Φ Mode",   true));

    const juce::StringArray noiseChoices { "White","Pink","Brown","Aureo","Quantum","Velvet" };
    p.push_back (std::make_unique<juce::AudioParameterChoice>(noiseMode,  "Noise Mode",  noiseChoices, 5));
    p.push_back (std::make_unique<juce::AudioParameterChoice>(noiseColor, "Noise Color", noiseChoices, 1));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(colorAmt,    "Color Amt", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.65f));

    // Aureo/Quantum advanced
    p.push_back (std::make_unique<juce::AudioParameterFloat>(aureoMix,   "Aureo Mix",   juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.5f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(aureoDecay, "Aureo Decay", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.3f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(aureoStride,"Aureo Stride",juce::NormalisableRange<float> (0.1f, 3.0f, 0.0001f, 0.5f), 1.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(aureoHarmo, "Aureo Harm",  juce::NormalisableRange<float> (1.0f, 24.0f, 1.0f), 12.0f));
    p.push_back (std::make_unique<juce::AudioParameterBool>(aureoVelvet, "Aureo Velvet", true));

    p.push_back (std::make_unique<juce::AudioParameterFloat>(qMix,     "Quantum Mix",   juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.55f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(qDetail,  "Quantum Detail",juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.7f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(qStride,  "Quantum Stride",juce::NormalisableRange<float> (0.1f, 3.0f, 0.0001f, 0.5f), 1.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(qBase,    "Quantum Base",  juce::NormalisableRange<float> (20.0f, 2000.0f, 0.01f, 0.4f), 220.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(qVelvet,  "Quantum Velvet",juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.12f));

    // VHS
    p.push_back (std::make_unique<juce::AudioParameterFloat>(vhsWow,     "VHS Wow",     juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.35f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(vhsFlutter, "VHS Flutter", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.25f));

    // Φ extended controls
    p.push_back (std::make_unique<juce::AudioParameterFloat>(phiRatio,    "Φ Ratio",    juce::NormalisableRange<float> (1.0f, 2.5f, 0.0001f, 0.45f), (float) aureo::kPhi));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(phiHeadB,    "Head b",     juce::NormalisableRange<float> (0.06f, 0.11f, 0.0001f), 0.0875f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(phiDistance, "Φ Distance", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.12f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(phiElev,     "Φ Elev",     juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.6f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(phiElevNotch,"Φ Elev Notch",juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.6f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(phiTorsoMix, "Torso Mix",  juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.12f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(phiTorsoMs,  "Torso ms",   juce::NormalisableRange<float> (0.1f, 4.0f, 0.0001f, 0.5f), 1.15f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(phiTorsoHp,  "Torso HP Hz",juce::NormalisableRange<float> (50.0f, 4000.0f, 0.01f, 0.4f), 500.0f));

    // pinna
    p.push_back (std::make_unique<juce::AudioParameterBool>(pinnaOn,    "Pinna",    false));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(pinnaDepth,"Pinna Depth", juce::NormalisableRange<float> (0.0f, 24.0f, 0.0001f), 12.0f));

    // spatial profile + CEFG
    p.push_back (std::make_unique<juce::AudioParameterChoice>(spatialProfile, "Spatial", juce::StringArray{ "Legacy", "UteroidΦ" }, 0));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(cefgGain, "CEFG Gain", juce::NormalisableRange<float> (0.0f, 4.0f, 0.0001f), 1.1f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(cefgMix,  "CEFG Mix",  juce::NormalisableRange<float> (0.0f, 4.0f, 0.0001f), 3.4f));

    // host sync
    p.push_back (std::make_unique<juce::AudioParameterBool>(syncEnable, "Sync", false));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(syncDivision, "Division", juce::NormalisableRange<float> (0.0625f, 32.0f, 0.0001f, 0.5f), 1.0f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(syncSlew,  "Sync Slew", juce::NormalisableRange<float> (0.0f, 5.0f, 0.0001f), 0.5f));

    // modal
    p.push_back (std::make_unique<juce::AudioParameterBool>(modalOn,    "Modal", false));
    p.push_back (std::make_unique<juce::AudioParameterBool>(modalMirror,"Modal Mirror", false));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(modalMix,  "Modal Mix", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.35f));
    p.push_back (std::make_unique<juce::AudioParameterFloat>(modalDecay,"Modal Decay", juce::NormalisableRange<float> (0.0f, 1.0f, 0.0001f), 0.5f));
    p.push_back (std::make_unique<juce::AudioParameterChoice>(modalPreset, "Modal Preset", juce::StringArray{ "Off","Wood","Metal","Glass" }, 1));

    return { p.begin(), p.end() };
}

// old inline Engine implementation removed (replaced by Mirror7Engine)
