#pragma once

#include <memory>
#include <ostream>
#include <string>

class Painter {
 public:
  Painter(std::ostream& out) : out_(out) {}
  virtual ~Painter() {}

  virtual void Prepare() const = 0;
  virtual void SetAxis(double xmin, double xmax, double ymin, double ymax) const = 0;
  virtual void DrawLine(double x1, double y1, double x2, double y2) const = 0;
  virtual void DrawDottedLine(double x1, double y1, double x2, double y2) const = 0;
  virtual void DrawArrow(double x1, double y1, double x2, double y2) const = 0;
  virtual void DrawPoint(double x, double y) const = 0;
  virtual void DrawText(double x, double y, const std::string& text) const = 0;
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
  virtual void SetAxis(double xmin, double xmax, double ymin, double ymax) const override;
  virtual void DrawLine(double x1, double y1, double x2, double y2) const override;
  virtual void DrawDottedLine(double x1, double y1, double x2, double y2) const override;
  virtual void DrawArrow(double x1, double y1, double x2, double y2) const override;
  virtual void DrawPoint(double x, double y) const override;
  virtual void DrawText(double x, double y, const std::string& text) const override;
  virtual void Finish() const override;
};

class AsymptotePainter : public Painter {
 public:
  AsymptotePainter(std::ostream& out) : Painter(out) {}
  virtual ~AsymptotePainter() {}

  virtual void Prepare() const override;
  virtual void SetAxis(double xmin, double xmax, double ymin, double ymax) const override;
  virtual void DrawLine(double x1, double y1, double x2, double y2) const override;
  virtual void DrawDottedLine(double x1, double y1, double x2, double y2) const override;
  virtual void DrawArrow(double x1, double y1, double x2, double y2) const override;
  virtual void DrawPoint(double x, double y) const override;
  virtual void DrawText(double x, double y, const std::string& text) const override;
  virtual void Finish() const override;
};

class PainterFactory {
 public:
  static std::unique_ptr<Painter> create(PainterBackend backend, std::ostream& out) {
    switch (backend) {
      case PainterBackend::kTikz:
        return std::make_unique<TikzPainter>(out);
      case PainterBackend::kAsymptote:
        return std::make_unique<AsymptotePainter>(out);
    }
    return nullptr;
  }
};