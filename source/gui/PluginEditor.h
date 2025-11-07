#pragma once

#include <juce_gui_extra/juce_gui_extra.h>
#include "plugin/PluginProcessor.h"
#include "config/ParamIDs.h"

class Mirror7AudioProcessorEditor : public juce::AudioProcessorEditor, private juce::Button::Listener
{
public:
    explicit Mirror7AudioProcessorEditor (Mirror7AudioProcessor& p);
    ~Mirror7AudioProcessorEditor() override = default;

    void paint (juce::Graphics&) override;
    void resized() override;
    void buttonClicked (juce::Button* b) override;

private:
    Mirror7AudioProcessor& processor;

    // Top bar
    juce::Label title { {}, "Mirror7" };
    juce::TextButton saveBtn { "Save Preset" }, loadBtn { "Load Preset" };
    juce::ComboBox oversample;
    std::unique_ptr<juce::AudioProcessorValueTreeState::ComboBoxAttachment> osA;

    // Scrollable content with all controls
    juce::Viewport viewport;
    std::unique_ptr<juce::Component> content;

    // Storage for attachments and widgets
    struct SliderItem { std::unique_ptr<juce::Label> label; std::unique_ptr<juce::Slider> slider; std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> att; };
    struct ToggleItem { std::unique_ptr<juce::ToggleButton> button; std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment> att; };
    struct ChoiceItem { std::unique_ptr<juce::Label> label; std::unique_ptr<juce::ComboBox> combo; std::unique_ptr<juce::AudioProcessorValueTreeState::ComboBoxAttachment> att; };

    std::vector<SliderItem> sliders;
    std::vector<ToggleItem> toggles;
    std::vector<ChoiceItem> choices;

    void addSlider (const juce::String& name, const juce::String& paramID);
    void addToggle (const juce::String& name, const juce::String& paramID);
    void addChoice (const juce::String& name, const juce::String& paramID, const juce::StringArray& items);
    void layoutContent();
};
