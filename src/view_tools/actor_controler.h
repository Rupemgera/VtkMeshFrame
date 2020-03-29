/**
 * @file actor_controler.h
 * @author Liao Yizhou (3130001008@zju.edu.cn)
 * @brief class for managing vtkAcotrs
 * @version 0.1
 * @date 2020-01-03
 *
 * @copyright Copyright (c) 2020
 *
 */

#pragma once
#include "view_defs.h"
#include "vtk_heads.h"

namespace viewtools {

// all subclass of ActorControler

class Actor3dControler {
public:
  enum class ActorType {
    BASE,
    MESH,
    SEGMENTS,
    POINTS,
    VECTORS,
    MOLECULES,
    LABEL
  };

  using ActorPtr = vtkSmartPointer<vtkActor>;

protected:
  std::string _name{""};

  ActorPtr _actor{nullptr};

  bool visibility_flag{true};

  ActorType _type{ActorType::BASE};

  // derived class should call this constructor and specify their real type
  Actor3dControler(std::string name, ActorPtr actor, ActorType type)
      : _name{name}, _actor{actor}, _type{type} {}

public:
  Actor3dControler(std::string name, ActorPtr actor);
  virtual ~Actor3dControler();

  std::string getName() const { return _name; }

  void setVisibility(bool visibility);

  bool getVisibility();

  virtual void setRenderSyle(int nRenderStyle);

  virtual void setOpacity(double opacity);

  virtual void setColor(Color color);

  virtual void setColor(double r, double g, double b);

  void setIntColor(int r, int g, int b);

  virtual std::tuple<int, int, int> getIntColor();

  virtual void setSize(double size);

  virtual double getSize();

  /**
   *@brief return the actual class type of the instance
   */
  ActorType getClassType();

  // tuple 1.visibility 2.opacity 3.size
  virtual std::tuple<bool, double, double> getStatus();

  ActorPtr get_actor() { return _actor; }

  virtual vtkProp *getProp() { return _actor.Get(); }

  vtkDataSet *get_data() { return _actor->GetMapper()->GetInput(); }
};

class MeshActorControler : public Actor3dControler {
public:
  struct {
    // whether edges are rendered
    bool edge_on{false};

    Color edge_color{0.0, 0.545, 0.0};

    // whether faces are rendered
    bool face_on{false};

    Color face_color{1.0, 1.0, 1.0};
  } render_status;

  MeshActorControler(std::string name, ActorPtr actor);

  /**
   * 1: edge 2:face 3:both
   */
  void setRenderSyle(int nRenderStyle);

  vtkPolyData *get_data() {
    vtkPolyDataMapper *mapper =
        static_cast<vtkPolyDataMapper *>(_actor->GetMapper());
    return mapper->GetInput();
  }
};

class PointsActorControler : public Actor3dControler {
public:
  struct {
    ////// properties of points //////////////

    float point_size{1.0};

    Color point_color{0.4157, 0.3529, 0.8039};
  } render_status;

  PointsActorControler(std::string name, ActorPtr actor);

  ////////////// function about points render //////////////////

  /**
   * set size of points
   */
  void setSize(double size) override;

  double getSize() override;

  std::tuple<bool, double, double> getStatus() override;
};

class SegmentActorControler : public Actor3dControler {
public:
  SegmentActorControler(std::string name, ActorPtr actor)
      : Actor3dControler(name, actor, ActorType::SEGMENTS) {}

  void setSize(double size) override;
};

class VectorActorControler : public Actor3dControler {
public:
  VectorActorControler(std::string name, ActorPtr actor)
      : Actor3dControler(name, actor, ActorType::VECTORS) {}

  // void setSize(double size) override;
};

class MoleculeActorControler : public Actor3dControler {
public:
  MoleculeActorControler(std::string name, ActorPtr actor)
      : Actor3dControler(name, actor, ActorType::MOLECULES) {}
};

class LabelActorControler : public Actor3dControler {
protected:
  vtkSmartPointer<vtkActor2D> _actor2d;

public:
  LabelActorControler(std::string name, vtkSmartPointer<vtkActor2D> actor)
      : Actor3dControler(name, nullptr, ActorType::LABEL), _actor2d{actor} {}
};

class Actor2dControler {
protected:
  std::string _name;

  vtkSmartPointer<vtkActor2D> _actor2d;

public:
  Actor2dControler(std::string name, vtkSmartPointer<vtkActor2D> actor)
      : _name{name}, _actor2d{actor} {}

  vtkProp *getProp() { return _actor2d.Get(); }

  std::string getName() { return _name; }

  void setVisibility(bool flag);

  bool getVisibility() { return _actor2d->GetVisibility(); }
};
} // namespace viewtools