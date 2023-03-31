#include <cassert>
#include <vector>

#include "Painter/Painter.h"

class AsymptoteNode {
 public:
  virtual ~AsymptoteNode() {}
  virtual std::string ToString() const = 0;
  friend std::ostream& operator<<(std::ostream& out, const AsymptoteNode& node);
};

std::ostream& operator<<(std::ostream& out, const AsymptoteNode& node) {
  out << node.ToString();
  return out;
}

class AsymptotePoint : public AsymptoteNode {
 public:
  AsymptotePoint(double x, double y) : x_(x), y_(y) {}
  virtual ~AsymptotePoint() = default;

  virtual std::string ToString() const override {
    return "(" + std::to_string(x_) + ", " + std::to_string(y_) + ")";
  }

 private:
  double x_;
  double y_;
};

class AsymptotePath : public AsymptoteNode {
 public:
  AsymptotePath(const std::vector<AsymptotePoint>& points) : points_(points) {}
  virtual ~AsymptotePath() = default;

  virtual std::string ToString() const override {
    assert(points_.size() >= 2);
    std::string path = points_[0].ToString();
    for (int i = 1; i < points_.size(); ++i) {
      path += " -- " + points_[i].ToString();
    }
    return path;
  }

 private:
  std::vector<AsymptotePoint> points_;
};

static std::string get_options(const std::vector<std::string>& options) {
  std::string options_str = "";

  if (!options.empty()) {
    options_str = ",";
    options_str += options[0];
    for (int i = 1; i < options.size(); ++i) {
      options_str += ", " + options[i];
    }
  }

  return options_str;
}

class AsymptoteDraw : public AsymptoteNode {
 public:
  AsymptoteDraw(const AsymptoteNode& node, const std::vector<std::string>& options)
      : node_(node), options_(options) {}
  virtual ~AsymptoteDraw() = default;

  virtual std::string ToString() const override {
    std::string options = get_options(options_);
    return "draw(" + node_.ToString() + options + ");\n";
  }

 private:
  const AsymptoteNode& node_;
  std::vector<std::string> options_;
};

class AsymptoteFilledDraw : public AsymptoteNode {
 public:
  AsymptoteFilledDraw(const AsymptoteNode& node, const std::vector<std::string>& options)
      : node_(node), options_(options) {}
  virtual ~AsymptoteFilledDraw() = default;

  virtual std::string ToString() const override {
    std::string options = get_options(options_);
    return "fill(" + node_.ToString() + options + ");\n";
  }

 private:
  const AsymptoteNode& node_;
  std::vector<std::string> options_;
};

class AsymptoteLine : public AsymptoteNode {
 public:
  AsymptoteLine(const std::vector<AsymptotePoint>& points, const std::vector<std::string>& options)
      : points_(points), options_(options) {}
  virtual ~AsymptoteLine() = default;

  virtual std::string ToString() const override {
    assert(points_.size() >= 2);
    AsymptotePath path(points_);
    AsymptoteDraw draw(path, options_);
    return draw.ToString();
  }

 private:
  std::vector<AsymptotePoint> points_;
  std::vector<std::string> options_;
};

class AsymptoteText : public AsymptoteNode {
 public:
  AsymptoteText(const AsymptoteNode& node, const std::string& text,
                const std::vector<std::string>& options)
      : node_(node), text_(text), options_(options) {}
  virtual ~AsymptoteText() = default;

  virtual std::string ToString() const override {
    std::string options = get_options(options_);
    return "label(\"" + text_ + "\", " + node_.ToString() + options + ");\n";
  }

 private:
  const AsymptoteNode& node_;
  std::string text_;
  std::vector<std::string> options_;
};

class AsymptoteCircle : public AsymptoteNode {
 public:
  AsymptoteCircle(const AsymptotePoint& center, const std::vector<std::string>& options)
      : center_(center), options_(options) {}
  virtual ~AsymptoteCircle() = default;

  virtual std::string ToString() const override {
    std::string options = get_options(options_);
    return "circle(" + center_.ToString() + options + ")";
  }

 private:
  AsymptotePoint center_;
  std::vector<std::string> options_;
};

void AsymptotePainter::Prepare() const {
  out_ << "import graph;\n";
  out_ << "size(500,500);\n";
};

void AsymptotePainter::SetAxis(double xmin, double xmax, double ymin, double ymax) const {
  out_ << AsymptoteLine({AsymptotePoint(xmin, 0), AsymptotePoint(xmax, 0)}, {});
  out_ << AsymptoteLine({AsymptotePoint(0, ymin), AsymptotePoint(0, ymax)}, {});
}

void AsymptotePainter::DrawLine(double x1, double y1, double x2, double y2) const {
  out_ << AsymptoteLine({AsymptotePoint(x1, y1), AsymptotePoint(x2, y2)}, {});
}

void AsymptotePainter::DrawArrow(double x1, double y1, double x2, double y2) const {
  out_ << AsymptoteLine({AsymptotePoint(x1, y1), AsymptotePoint(x2, y2)}, {"red", "Arrow"});
}

void AsymptotePainter::DrawDottedLine(double x1, double y1, double x2, double y2) const {
  out_ << AsymptoteLine({AsymptotePoint(x1, y1), AsymptotePoint(x2, y2)}, {"dotted"});
}

void AsymptotePainter::DrawPoint(double x, double y) const {
  out_ << AsymptoteFilledDraw(AsymptoteCircle(AsymptotePoint(x, y), {"0.1"}), {});
}

void AsymptotePainter::DrawText(double x, double y, const std::string& text) const {
  out_ << AsymptoteText(AsymptotePoint(x, y), text, {"fontsize(20pt)"});
}

void AsymptotePainter::Finish() const {
  return;  // nothing to do
}
