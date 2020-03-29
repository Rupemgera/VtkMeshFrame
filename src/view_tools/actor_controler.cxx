#include "actor_controler.h"

namespace viewtools { /************************* ActorControler begin
                       *************************/

Actor3dControler::Actor3dControler(std::string name, ActorPtr actor)
    : _name(name), _actor(actor) {}

Actor3dControler::~Actor3dControler() {}

void Actor3dControler::setVisibility(bool visibility) {
  if (_actor == nullptr)
    return;
  visibility_flag = visibility;
  _actor->SetVisibility(visibility);
}

/**
        surface type : points(0), wireframe(1) or surface(2)
*/
void Actor3dControler::setRenderSyle(int nRenderStyle) {
  if (_actor == nullptr)
    return;
  _actor->GetProperty()->SetRepresentation(nRenderStyle);
}

void Actor3dControler::setColor(Color color) {
  if (_actor == nullptr)
    return;
  _actor->GetProperty()->SetColor(color.data());
}

std::tuple<int, int, int> Actor3dControler::getIntColor() {
  double *c = _actor->GetProperty()->GetColor();
  int r = static_cast<int>(c[0] * 255.0);
  int g = static_cast<int>(c[1] * 255.0);
  int b = static_cast<int>(c[2] * 255.0);
  return std::tuple<int, int, int>(r, g, b);
}

void Actor3dControler::setColor(double r, double g, double b) {
  if (_actor == nullptr)
    return;
  _actor->GetProperty()->SetColor(r, g, b);
}

void Actor3dControler::setIntColor(int r, int g, int b) {
  this->setColor(r / 255.0, g / 255.0, b / 255.0);
}

void Actor3dControler::setSize(double size) {
  if (_actor == nullptr)
    return;
  _actor->GetProperty()->SetLineWidth(size);
}

double Actor3dControler::getSize() {
  if (_actor == nullptr)
    return 0.0;
  return _actor->GetProperty()->GetLineWidth();
}

bool Actor3dControler::getVisibility() { return visibility_flag; }

Actor3dControler::ActorType Actor3dControler::getClassType() { return _type; }

std::tuple<bool, double, double> Actor3dControler::getStatus() {
  double opacity = _actor->GetProperty()->GetOpacity();
  double size = _actor->GetProperty()->GetLineWidth();
  return std::tuple<bool, double, double>(visibility_flag, opacity, size);
}

void Actor3dControler::setOpacity(double opacity) {
  if (_actor == nullptr)
    return;
  _actor->GetProperty()->SetOpacity(opacity);
}

/************************* MeshActorControler  begein *************************/

MeshActorControler::MeshActorControler(std::string name, ActorPtr actor)
    : Actor3dControler(name, actor, ActorType::MESH) {
  _actor->GetProperty()->SetColor(render_status.face_color.data());
  _actor->GetProperty()->SetEdgeColor(render_status.edge_color.data());
}

void MeshActorControler::setRenderSyle(int nRenderStyle) {
  if (nRenderStyle & 2) {
    _actor->VisibilityOn();
    // surface type : points(0), wireframe(1) or surface(2)
    _actor->GetProperty()->SetRepresentation(2);
    render_status.face_on = true;
    _actor->GetProperty()->SetColor(render_status.face_color.data());
    // also render edges
    if (nRenderStyle & 1) {
      _actor->GetProperty()->SetEdgeVisibility(true);
      _actor->GetProperty()->SetEdgeColor(render_status.edge_color.data());
      render_status.edge_on = true;
    } else {
      _actor->GetProperty()->SetEdgeVisibility(false);
      render_status.edge_on = false;
    }
  } else if (nRenderStyle & 1) // only render edges
  {
    _actor->VisibilityOn();
    _actor->GetProperty()->SetRepresentationToWireframe();
    _actor->GetProperty()->SetColor(render_status.edge_color.data());
    render_status.face_on = false;
    render_status.edge_on = true;
  } else {
    _actor->VisibilityOff();
    render_status.face_on = false;
    render_status.edge_on = false;
  }
}

/************************* MeshActorControler  begein *************************/

PointsActorControler::PointsActorControler(std::string name, ActorPtr actor)
    : Actor3dControler(name, actor, ActorType::POINTS) {
  _actor->GetProperty()->SetColor(render_status.point_color.data());
}

void PointsActorControler::setSize(double size) {
  render_status.point_size = size;
  _actor->GetProperty()->SetPointSize(size);
}

double PointsActorControler::getSize() {
  return _actor->GetProperty()->GetPointSize();
}

std::tuple<bool, double, double> PointsActorControler::getStatus() {
  double opacity = _actor->GetProperty()->GetOpacity();
  double size = _actor->GetProperty()->GetPointSize();
  return std::tuple<bool, double, double>(visibility_flag, opacity, size);
}

/************************* SegmentActorControler  *************************/

void SegmentActorControler::setSize(double size) {
  _actor->GetProperty()->SetLineWidth(size);
}

void Actor2dControler::setVisibility(bool flag) {
  if (_actor2d != nullptr) {
    _actor2d->SetVisibility(flag);
  }
}

} // namespace viewtools