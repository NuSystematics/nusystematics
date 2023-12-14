#pragma once

#include "nusystematics/responsecalculators/TemplateResponseCalculatorBase.hh"

namespace nusyst {

class TemplateResponseQ0Q3EAvail : public TemplateResponse3DDiscrete {
public:
  std::string GetCalculatorName() const { return "TemplateResponseQ0Q3EAvail"; }
};

typedef TemplateResponseQ0Q3EAvail FSILikeEAvailSmearing_ReWeight;
} // namespace nusyst
