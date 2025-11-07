#include "gui/PluginEditor.h"
#include "config/ParamIDs.h"
using namespace ParamIDs;

static void configKnob (juce::Slider& s, juce::Label& l, const juce::String& name)
{
    l.setText (name, juce::dontSendNotification);
    l.setJustificationType (juce::Justification::centredLeft);
    s.setSliderStyle (juce::Slider::RotaryHorizontalVerticalDrag);
    s.setTextBoxStyle (juce::Slider::TextBoxBelow, false, 60, 18);
}

Mirror7AudioProcessorEditor::Mirror7AudioProcessorEditor (Mirror7AudioProcessor& p)
    : juce::AudioProcessorEditor (&p), processor (p)
{
    setSize (760, 540);

    addAndMakeVisible (title);
    title.setJustificationType (juce::Justification::centred);
    title.setFont (juce::Font (juce::FontOptions (18.0f, juce::Font::bold)));

    auto& vts = processor.getVTS();

    // Top bar buttons
    addAndMakeVisible (saveBtn); addAndMakeVisible (loadBtn);
    saveBtn.addListener (this); loadBtn.addListener (this);

    addAndMakeVisible (oversample);
    oversample.addItem ("1x", 1); oversample.addItem ("2x", 2); oversample.addItem ("4x", 3);
    osA = std::make_unique<juce::AudioProcessorValueTreeState::ComboBoxAttachment> (vts, osFactor, oversample);

    // Content container in a viewport (scrollable)
    content = std::make_unique<juce::Component>();
    addAndMakeVisible (viewport);
    viewport.setViewedComponent (content.get(), false);

    // Add controls (grouped logically)
    addSlider ("Output", outGain);
    addChoice ("Noise Mode", noiseMode, juce::StringArray());
    addChoice ("Noise Color", noiseColor, juce::StringArray());
    addSlider ("Color Amt", colorAmt);
    addChoice ("Oversample", osFactor, juce::StringArray()); // already shown; keeps state in content

    // Stereo / pan
    addSlider ("Pan", pan); addSlider ("Width", width);
    addToggle ("Autopan", autopan); addSlider ("Autopan Hz", autopanSp);

    // ITD/ILD + spatial guard
    addSlider ("ITD (us)", itdUs); addSlider ("ILD (dB)", ildDb);
    addSlider ("Min Deg", spatMinDeg); addSlider ("Min ms", spatMinMs);
    addSlider ("IPD", spatIPD); addSlider ("Shadow", spatShadow); addToggle ("Mirror", spatMirror);

    // Grain timing
    addSlider ("Rate (Hz)", rateHz); addSlider ("Base Len (ms)", baseMs); addSlider ("Len Φ", lenPhi); addSlider ("Hemis Coupling", hemisCoupling);

    // Envelope
    addSlider ("Env A", envA); addSlider ("Env D", envD); addSlider ("Env S", envS); addSlider ("Env R", envR);

    // Dialogue
    addToggle ("Dialogue", dialogOn); addSlider ("Dial Str", dialogStr); addSlider ("Dial Mem", dialogMem); addSlider ("Dial ΦMix", dialogPhi);

    // Glitch/Crush
    addSlider ("Glitch", glitchMix); addSlider ("SR Crush", srCrush); addSlider ("BitCrush", bitCrush);
    addSlider ("VHS Wow", vhsWow); addSlider ("VHS Flutter", vhsFlutter);

    // Noise advanced
    addSlider ("Velvet", velvetAmt);
    addSlider ("Aureo Mix", aureoMix); addSlider ("Aureo Decay", aureoDecay); addSlider ("Aureo Stride", aureoStride); addSlider ("Aureo Harm", aureoHarmo); addToggle ("Aureo Velvet", aureoVelvet);
    addSlider ("Quantum Mix", qMix); addSlider ("Quantum Detail", qDetail); addSlider ("Quantum Stride", qStride); addSlider ("Quantum Base", qBase); addSlider ("Quantum Velvet", qVelvet);

    // Φ controls
    addToggle ("Φ Mode", phiMode); addToggle ("Φ Pan", phiPan);
    addSlider ("Φ Ratio", phiRatio); addSlider ("Head b", phiHeadB); addSlider ("Φ Distance", phiDistance);
    addSlider ("Φ Elev", phiElev); addSlider ("Φ Elev Notch", phiElevNotch); addSlider ("Torso Mix", phiTorsoMix); addSlider ("Torso ms", phiTorsoMs); addSlider ("Torso HP Hz", phiTorsoHp);

    // Modal
    addToggle ("Modal", modalOn); addToggle ("Modal Mirror", modalMirror); addSlider ("Modal Mix", modalMix); addSlider ("Modal Decay", modalDecay);
    addChoice ("Modal Preset", modalPreset, juce::StringArray());

    // Pinna
    addToggle ("Pinna", pinnaOn); addSlider ("Pinna Depth", pinnaDepth);

    // Spatial v2
    addChoice ("Spatial", spatialProfile, juce::StringArray()); addSlider ("CEFG Gain", cefgGain); addSlider ("CEFG Mix", cefgMix);

    // Sync
    addToggle ("Sync", syncEnable); addSlider ("Division", syncDivision); addSlider ("Sync Slew", syncSlew);

    layoutContent();
}

void Mirror7AudioProcessorEditor::paint (juce::Graphics& g)
{
    g.fillAll (juce::Colour::fromRGB (18, 18, 26));
    g.setColour (juce::Colours::white);
    g.drawRect (getLocalBounds());
}

void Mirror7AudioProcessorEditor::resized()
{
    auto r = getLocalBounds().reduced (12);
    auto top = r.removeFromTop (36);
    title.setBounds (top.removeFromLeft (160));
    saveBtn.setBounds (top.removeFromLeft (140).reduced (4));
    loadBtn.setBounds (top.removeFromLeft (140).reduced (4));
    oversample.setBounds (top.removeFromRight (120).reduced (8, 4));

    viewport.setBounds (r);
    layoutContent();
}

void Mirror7AudioProcessorEditor::buttonClicked (juce::Button* b)
{
    if (b == &saveBtn)
    {
        auto chooser = std::make_shared<juce::FileChooser> ("Save Mirror7 Preset", juce::File(), "*.mir7preset");
        juce::Component::SafePointer<Mirror7AudioProcessorEditor> safeThis (this);
        chooser->launchAsync (juce::FileBrowserComponent::saveMode | juce::FileBrowserComponent::canSelectFiles,
                              [safeThis, chooser](const juce::FileChooser& fc)
                              {
                                  if (auto* editor = safeThis.getComponent())
                                  {
                                      auto f = fc.getResult();
                                      if (f == juce::File()) return;
                                      if (! f.hasFileExtension ("mir7preset")) f = f.withFileExtension ("mir7preset");
                                      juce::MemoryOutputStream mos;
                                      editor->processor.getVTS().state.writeToStream (mos);
                                      f.replaceWithData (mos.getData(), mos.getDataSize());
                                  }
                              });
    }
    else if (b == &loadBtn)
    {
        auto chooser = std::make_shared<juce::FileChooser> ("Load Mirror7 Preset", juce::File(), "*.mir7preset");
        juce::Component::SafePointer<Mirror7AudioProcessorEditor> safeThis (this);
        chooser->launchAsync (juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles,
                              [safeThis, chooser](const juce::FileChooser& fc)
                              {
                                  if (auto* editor = safeThis.getComponent())
                                  {
                                      auto f = fc.getResult();
                                      if (f == juce::File()) return;
                                      juce::FileInputStream fis (f);
                                      if (fis.openedOk())
                                      {
                                          auto vt = juce::ValueTree::readFromStream (fis);
                                          if (vt.isValid())
                                              editor->processor.getVTS().replaceState (vt);
                                      }
                                  }
                              });
    }
}

void Mirror7AudioProcessorEditor::addSlider (const juce::String& name, const juce::String& paramID)
{
    auto& vts = processor.getVTS();
    SliderItem it;
    it.label = std::make_unique<juce::Label>();
    it.slider = std::make_unique<juce::Slider>();
    configKnob (*it.slider, *it.label, name);
    content->addAndMakeVisible (*it.label);
    content->addAndMakeVisible (*it.slider);
    it.att = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment> (vts, paramID, *it.slider);
    sliders.push_back (std::move (it));
}

void Mirror7AudioProcessorEditor::addToggle (const juce::String& name, const juce::String& paramID)
{
    auto& vts = processor.getVTS();
    ToggleItem it;
    it.button = std::make_unique<juce::ToggleButton> (name);
    content->addAndMakeVisible (*it.button);
    it.att = std::make_unique<juce::AudioProcessorValueTreeState::ButtonAttachment> (vts, paramID, *it.button);
    toggles.push_back (std::move (it));
}

void Mirror7AudioProcessorEditor::addChoice (const juce::String& name, const juce::String& paramID, const juce::StringArray&)
{
    auto& vts = processor.getVTS();
    ChoiceItem it;
    it.label = std::make_unique<juce::Label>();
    it.combo = std::make_unique<juce::ComboBox>();
    it.label->setText (name, juce::dontSendNotification);
    it.label->setJustificationType (juce::Justification::centredLeft);
    content->addAndMakeVisible (*it.label);
    content->addAndMakeVisible (*it.combo);
    if (auto* param = vts.getParameter (paramID))
    {
        if (auto* apc = dynamic_cast<juce::AudioParameterChoice*> (param))
        {
            const auto& opts = apc->choices;
            juce::StringArray items;
            for (int i = 0; i < opts.size(); ++i) items.add (opts[i]);
            for (int i = 0; i < items.size(); ++i) it.combo->addItem (items[i], i + 1);
        }
    }
    it.att = std::make_unique<juce::AudioProcessorValueTreeState::ComboBoxAttachment> (vts, paramID, *it.combo);
    choices.push_back (std::move (it));
}

void Mirror7AudioProcessorEditor::layoutContent()
{
    if (content == nullptr) return;
    auto area = viewport.getLocalBounds().withWidth (getWidth() - 24); // leave margins
    const int cols = 3;
    const int cellW = area.getWidth() / cols;
    const int controlH = 68;
    int x = 0, y = 0, col = 0;

    auto placeSlider = [&] (SliderItem& s)
    {
        auto r = juce::Rectangle<int> (col * cellW + 8, y + 8, cellW - 16, controlH - 16);
        s.slider->setBounds (r);
        s.label->setBounds (r.withHeight (18).withY (y));
        if (++col >= cols) { col = 0; y += controlH; }
    };
    auto placeToggle = [&] (ToggleItem& t)
    {
        t.button->setBounds (col * cellW + 12, y + 20, cellW - 24, 28);
        if (++col >= cols) { col = 0; y += controlH; }
    };
    auto placeChoice = [&] (ChoiceItem& c)
    {
        c.label->setBounds (col * cellW + 10, y + 4, cellW - 20, 18);
        c.combo->setBounds (col * cellW + 10, y + 24, cellW - 20, 24);
        if (++col >= cols) { col = 0; y += controlH; }
    };

    for (auto& s : sliders) placeSlider (s);
    for (auto& t : toggles) placeToggle (t);
    for (auto& c : choices) placeChoice (c);

    const int totalH = y + controlH + 12;
    content->setBounds (8, 0, area.getWidth(), totalH);
}
