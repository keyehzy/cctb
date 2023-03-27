#include <cassert>
#include <vector>

#include "Painter/Painter.h"

class TikzNode {
 public:
  virtual ~TikzNode() = default;
  virtual std::string ToString() const = 0;
  friend std::ostream& operator<<(std::ostream& out, const TikzNode& node);
};

std::ostream& operator<<(std::ostream& out, const TikzNode& node) {
  out << node.ToString();
  return out;
}

class TikzPoint : public TikzNode {
 public:
  TikzPoint(float x, float y) : x_(x), y_(y) {}
  virtual ~TikzPoint() = default;

  virtual std::string ToString() const override {
    return "(" + std::to_string(x_) + "," + std::to_string(y_) + ")";
  }

 private:
  float x_;
  float y_;
};

class TikzPath : public TikzNode {
 public:
  TikzPath(const std::vector<TikzPoint>& points) : points_(points) {}
  virtual ~TikzPath() = default;

  virtual std::string ToString() const override {
    assert(points_.size() >= 2);
    std::string path = points_[0].ToString();
    for (size_t i = 1; i < points_.size(); ++i) {
      path += " -- " + points_[i].ToString();
    }
    return path;
  }

 private:
  std::vector<TikzPoint> points_;
};

std::string get_options(const std::vector<std::string>& options) {
  std::string options_str = "";

  if (!options.empty()) {
    options_str += "[" + options[0];
    for (size_t i = 1; i < options.size(); ++i) {
      options_str += "," + options[i];
    }
    options_str += "]";
  }

  return options_str;
}

class TikzDraw : public TikzNode {
 public:
  TikzDraw(const TikzNode& node, std::vector<std::string> options)
      : node_(node), options_(options) {}
  virtual ~TikzDraw() = default;

  virtual std::string ToString() const override {
    std::string options = get_options(options_);
    return "\\draw" + options + " " + node_.ToString() + ";\n";
  }

 private:
  const TikzNode& node_;
  std::vector<std::string> options_;
};

class TikzFilledDraw : public TikzNode {
 public:
  TikzFilledDraw(const TikzNode& node, const TikzNode& obj, const std::vector<std::string>& options)
      : node_(node), obj_(obj), options_(options) {}
  virtual ~TikzFilledDraw() = default;

  virtual std::string ToString() const override {
    std::string options = get_options(options_);
    return "\\fill" + options + " " + node_.ToString() + obj_.ToString() + ";\n";
  }

 private:
  const TikzNode& node_;
  const TikzNode& obj_;
  std::vector<std::string> options_;
};

class TikzLine : public TikzNode {
 public:
  TikzLine(const std::vector<TikzPoint>& points, const std::vector<std::string>& options)
      : children_(points), options_(options) {}
  virtual ~TikzLine() = default;

  virtual std::string ToString() const override {
    assert(children_.size() >= 2);
    TikzPath path(children_);
    TikzDraw draw(path, options_);
    return draw.ToString();
  }

 private:
  std::vector<TikzPoint> children_;
  std::vector<std::string> options_;
};

class TikzText : public TikzNode {
 public:
  TikzText(const TikzNode& node, const std::string& text, std::vector<std::string> options)
      : node_(node), text_(text), options_(options) {}
  virtual ~TikzText() = default;

  virtual std::string ToString() const override {
    std::string options = "";

    if (!options_.empty()) {
      options += "[" + options_[0];
      for (size_t i = 1; i < options_.size(); ++i) {
        options += "," + options_[i];
      }
      options += "]";
    }

    return "\\node at " + node_.ToString() + " " + options + " {" + text_ + "};\n";
  }

 private:
  const TikzNode& node_;
  std::string text_;
  std::vector<std::string> options_;
};

class TikzCircle : public TikzNode {
 public:
  TikzCircle(float radius) : radius_(radius) {}
  virtual ~TikzCircle() = default;

  virtual std::string ToString() const override {
    return " circle (" + std::to_string(radius_) + ")";
  }

 private:
  float radius_;
};

void TikzPainter::Prepare() const {
  out_ << "\\documentclass[tikz]{standalone}\n";
  out_ << "\\begin{document}\n";
  out_ << "\\begin{tikzpicture}\n";
}

void TikzPainter::SetAxis(float xmin, float xmax, float ymin, float ymax) const {
  out_ << TikzLine({TikzPoint(xmin, 0), TikzPoint(xmax, 0)}, {});
  out_ << TikzLine({TikzPoint(0, ymin), TikzPoint(0, ymax)}, {});
}

void TikzPainter::DrawLine(float x1, float y1, float x2, float y2) const {
  out_ << TikzLine({TikzPoint(x1, y1), TikzPoint(x2, y2)}, {});
}

void TikzPainter::DrawArrow(float x1, float y1, float x2, float y2) const {
  out_ << TikzLine({TikzPoint(x1, y1), TikzPoint(x2, y2)}, {"->", "red"});
}

void TikzPainter::DrawDottedLine(float x1, float y1, float x2, float y2) const {
  out_ << TikzLine({TikzPoint(x1, y1), TikzPoint(x2, y2)}, {"dotted"});
}

void TikzPainter::DrawPoint(float x, float y) const {
  out_ << TikzFilledDraw(TikzPoint(x, y), TikzCircle(0.1), {"black"});
}

void TikzPainter::DrawText(float x, float y, const std::string& text) const {
  out_ << TikzText(TikzPoint(x, y), text, {"font=\\small"});
}

void TikzPainter::Finish() const {
  out_ << "\\end{tikzpicture}\n";
  out_ << "\\end{document}\n";
}