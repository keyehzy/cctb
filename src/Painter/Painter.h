#pragma once

#include <memory>
#include <ostream>
#include <string>

class Painter {
 public:
  Painter(std::ostream& out) : out_(out) {}
  virtual ~Painter() {}

  virtual void Prepare() const = 0;
  virtual void DrawLine(float x1, float y1, float x2, float y2) const = 0;
  virtual void DrawArrow(float x1, float y1, float x2, float y2) const = 0;
  virtual void DrawPoint(float x, float y) const = 0;
  virtual void DrawText(float x, float y, const std::string& text) const = 0;
  virtual void Finish() const = 0;

 protected:
  std::ostream& out_;
};

enum class PainterBackend {
  kTikz,
  kAsymptote,
};

class TikzPainter : public Painter {
 public:
  TikzPainter(std::ostream& out) : Painter(out) {}
  virtual ~TikzPainter() {}

  virtual void Prepare() const override;
  virtual void DrawLine(float x1, float y1, float x2, float y2) const override;
  virtual void DrawArrow(float x1, float y1, float x2, float y2) const override;
  virtual void DrawPoint(float x, float y) const override;
  virtual void DrawText(float x, float y,
                        const std::string& text) const override;
  virtual void Finish() const override;
};

class AsymptotePainter : public Painter {
 public:
  AsymptotePainter(std::ostream& out) : Painter(out) {}
  virtual ~AsymptotePainter() {}

  virtual void Prepare() const override;
  virtual void DrawLine(float x1, float y1, float x2, float y2) const override;
  virtual void DrawArrow(float x1, float y1, float x2, float y2) const override;
  virtual void DrawPoint(float x, float y) const override;
  virtual void DrawText(float x, float y,
                        const std::string& text) const override;
  virtual void Finish() const override;
};

class PainterFactory {
 public:
  static std::unique_ptr<Painter> create(PainterBackend backend,
                                         std::ostream& out) {
    switch (backend) {
      case PainterBackend::kTikz:
        return std::make_unique<TikzPainter>(out);
      case PainterBackend::kAsymptote:
        return std::make_unique<AsymptotePainter>(out);
    }
    return nullptr;
  }
};