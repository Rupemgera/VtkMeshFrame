#include "vtk_wrapper.h"

namespace viewtools {
/********************************/

vtkSmartPointer<vtkLookupTable>
get_lookuptable(double min, double max)
{ // Create a lookup table to share between the mapper and the
  // scalarbar
  vtkNew<vtkLookupTable> hueLut;
  hueLut->SetTableRange(min, max);
  // hueLut->SetHueRange(0.6, 1);
  // hueLut->SetSaturationRange(0.8, 1);
  // hueLut->SetValueRange(0.9, 0.9);
  hueLut->Build();
  return hueLut;
}

/********************************/

vtkSmartPointer<vtkPolyData>
VtkWrapper::processPolyData(const RowMatrix<double>& points,
                            const RowMatrix<IndexType>& polys)
{
  /* insert vertices */
  vtkNew<vtkPoints> nodes;
  size_t n_vertices = points.rows();
  nodes->GetData()->Allocate(n_vertices);
  for (int i = 0; i < n_vertices; ++i) {
    nodes->InsertPoint(i, points.row(i).data());
  }

  /* insert polys */
  vtkNew<vtkCellArray> cells;
  size_t cell_n = polys.cols();
  size_t n_faces = polys.rows();
  //每个单元有1+cell_n个数据 1个存储顶点个数cell_n，cell_n个存储面片顶点的标号
  cells->GetData()->Allocate((1 + cell_n) * n_faces);
  for (auto i = 0; i < n_faces; ++i) {
    cells->InsertNextCell(cell_n, polys.row(i).data());
  }

  /* form poly data */

  vtkNew<vtkPolyData> data;
  data->SetPoints(nodes);
  if (cell_n == 1)
    data->SetVerts(cells);
  else if (cell_n == 2)
    data->SetLines(cells);
  else
    data->SetPolys(cells);
  return data;
}

vtkSmartPointer<vtkPolyData>
VtkWrapper::processPointData(const RowMatrix<double>& points)
{
  /* insert vertices */
  size_t n_vertices = points.rows();
  auto nodes = convert_points(points);

  /* insert polys */
  vtkNew<vtkCellArray> cells;
  //每个单元有1+cell_n个数据 1个存储顶点个数cell_n，cell_n个存储面片顶点的标号
  cells->GetData()->Allocate(n_vertices);
  for (size_t i = 0; i < n_vertices; ++i) {
    cells->InsertNextCell(1);
    cells->InsertCellPoint(i);
  }

  /* form poly data */

  vtkNew<vtkPolyData> data;
  data->SetPoints(nodes);
  data->SetVerts(cells);
  return data;
}

vtkSmartPointer<vtkPoints>
VtkWrapper::convert_points(const RowMatrix<double>& points)
{
  assert(points.cols() == 3);
  vtkNew<vtkPoints> nodes;
  size_t n_vertices = points.rows();
  nodes->GetData()->Allocate(n_vertices);
  for (int i = 0; i < n_vertices; ++i) {
    nodes->InsertPoint(i, points.row(i).data());
  }
  return nodes;
}

vtkSmartPointer<vtkPoints>
VtkWrapper::convert_points(const std::vector<Eigen::Vector3d>& points)
{
  vtkNew<vtkPoints> nodes;
  size_t n_vertices = points.size();
  nodes->GetData()->Allocate(n_vertices);
  for (int i = 0; i < n_vertices; ++i) {
    nodes->InsertPoint(i, points[i].data());
  }
  return nodes;
}

vtkSmartPointer<vtkDoubleArray>
VtkWrapper::convert_arrary(const RowMatrix<double>& matrix)
{
  vtkNew<vtkDoubleArray> array;
  size_t tuple_size = matrix.cols();
  size_t n = matrix.rows();
  array->SetNumberOfComponents(static_cast<int>(tuple_size));
  array->SetNumberOfTuples(n);
  for (size_t i = 0; i < n; i++) {
    array->InsertTuple(i, matrix.row(i).data());
  }
  return array;
}

vtkSmartPointer<vtkDoubleArray>
VtkWrapper::convert_arrary(const std::vector<double>& scalars)
{
  vtkNew<vtkDoubleArray> array;
  size_t n = scalars.size();
  array->Allocate(n);
  for (size_t i = 0; i < n; i++) {
    array->InsertNextValue(scalars[i]);
  }
  return array;
}

vtkSmartPointer<vtkDoubleArray>
VtkWrapper::convert_vectors(const std::vector<Eigen::Vector3d>& vecs)
{
  vtkNew<vtkDoubleArray> array;
  size_t n = vecs.size();
  array->SetNumberOfComponents(3);
  array->SetNumberOfTuples(n);
  for (size_t i = 0; i < n; i++) {
    array->InsertTuple(i, vecs[i].data());
  }
  return array;
}

vtkSmartPointer<vtkDoubleArray>
VtkWrapper::convert_vectors(const RowMatrix<double>& vecs)
{
  vtkNew<vtkDoubleArray> array;
  size_t n = vecs.rows();
  array->SetNumberOfComponents(vecs.cols());
  array->SetNumberOfTuples(n);
  for (size_t i = 0; i < n; i++) {
    array->InsertTuple(i, vecs.row(i).data());
  }
  return array;
}

std::vector<vtkSmartPointer<vtkDoubleArray>>
VtkWrapper::convert_frames(const std::vector<Eigen::Matrix<double, 3, 3>>& vecs)
{
  std::vector<vtkSmartPointer<vtkDoubleArray>> res;
  res.clear();
  for (int dim = 0; dim < 3; ++dim) {
    vtkNew<vtkDoubleArray> array;
    size_t n = vecs.size();
    array->SetNumberOfComponents(3);
    array->SetNumberOfTuples(n);
    for (size_t i = 0; i < n; i++) {
      array->InsertTuple(i, vecs[i].col(dim).data());
    }
    res.push_back(array);
  }
  return res;
}

void
VtkWrapper::init_preset_colors()
{
  _preset_colors.clear();
  _preset_colors.push_back({ 1.0, 0.0, 0.0 });
  _preset_colors.push_back({ 0.0, 1.0, 0.0 });
  _preset_colors.push_back({ 0.0, 0.0, 1.0 });
  _preset_colors.push_back({ 1.0, 1.0, 0.0 });
  _preset_colors.push_back({ 1.0, 0.0, 1.0 });
  _preset_colors.push_back({ 0.0, 1.0, 1.0 });
  _preset_colors.push_back({ 0.0, 0.0, 0.0 });
}

/******************** public *******************************/

bool
VtkWrapper::drawArrows(std::string name,
                       const RowMatrix<double>& points,
                       const RowMatrix<double>& vectors)
{
  auto locs = convert_points(points);
  auto dirs = convert_vectors(vectors);

  vtkNew<vtkPolyData> pd;
  pd->SetPoints(locs);
  pd->GetPointData()->SetVectors(dirs);

  vtkNew<vtkArrowSource> arrow;
  vtkNew<vtkGlyph3D> glyph;
  glyph->SetSourceConnection(arrow->GetOutputPort());
  glyph->SetInputData(pd);
  // glyph->SetScaleFactor(1.0);
  glyph->SetScaleModeToScaleByVector();
  glyph->SetVectorModeToUseVector();

  // Mapper
  vtkNew<vtkPolyDataMapper> mapper;
  // mapper->SetInputConnection(input);
  mapper->SetInputConnection(glyph->GetOutputPort());

  // Actor in scene _actor->SetMapper(mapper)
  vtkNew<vtkActor> actor;

  // Actor in scene
  actor->SetMapper(mapper);

  auto ac = new VectorActorControler(name, actor);
  ac->setColor({ 1.0, 0.0, 1.0 });

  // actor->GetProperty()->SetColor(color.data());

  insert(ac);
  return true;
}

bool
VtkWrapper::drawArrows(std::string name,
                       const std::vector<Eigen::Vector3d>& points,
                       const std::vector<Eigen::Vector3d>& vectors)
{
  auto locs = convert_points(points);
  auto dirs = convert_vectors(vectors);

  vtkNew<vtkPolyData> pd;
  pd->SetPoints(locs);
  pd->GetPointData()->SetVectors(dirs);

  vtkNew<vtkArrowSource> arrow;
  vtkNew<vtkGlyph3D> glyph;
  glyph->SetSourceConnection(arrow->GetOutputPort());
  glyph->SetInputData(pd);
  // glyph->SetScaleFactor(1.0);
  glyph->SetScaleModeToScaleByVector();
  glyph->SetVectorModeToUseVector();

  // Mapper
  vtkNew<vtkPolyDataMapper> mapper;
  // mapper->SetInputConnection(input);
  mapper->SetInputConnection(glyph->GetOutputPort());

  // Actor in scene _actor->SetMapper(mapper)
  vtkNew<vtkActor> actor;

  // Actor in scene
  actor->SetMapper(mapper);

  auto ac = new VectorActorControler(name, actor);
  ac->setColor({ 1.0, 0.0, 1.0 });

  // actor->GetProperty()->SetColor(color.data());

  insert(ac);
  return true;
}

bool
VtkWrapper::drawFaces(std::string name,
                      const RowMatrix<double>& points,
                      const RowMatrix<IndexType>& faces)
{
  auto data = processPolyData(points, faces);

  /***** mapper *****/

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputData(data);

  /********** actors *************/

  vtkNew<vtkActor> face_actor;
  face_actor->SetMapper(mapper);

  auto ac = new MeshActorControler(name, face_actor);
  insert(ac);
  return true;
}

bool
VtkWrapper::drawFrames(std::string name,
                       const RowMatrix<double>& points,
                       const std::vector<Eigen::Matrix<double, 3, 3>>& frames,
                       double scale_factor,
                       double line_width)

{
  auto locs = convert_points(points);
  auto dirs = convert_frames(frames);
  Color colour[3] = { { 255, 0, 0 }, { 0, 255, 0 }, { 0, 0, 255 } };
  std::string label[3] = { "_major", "_middle", "_minor" };

  for (int dim = 0; dim < 3; ++dim) {
    vtkNew<vtkPolyData> ps;
    ps->SetPoints(locs);
    ps->GetPointData()->SetVectors(dirs[dim]);

    vtkNew<vtkHedgeHog> hedgehog;
    hedgehog->SetInputData(ps);
    hedgehog->SetScaleFactor(scale_factor);

    // Mapper
    vtkNew<vtkPolyDataMapper> mapper;
    // mapper->SetInputConnection(input);
    mapper->SetInputConnection(hedgehog->GetOutputPort());

    // Actor in scene _actor->SetMapper(mapper)
    vtkNew<vtkActor> actor;
    // Actor in scene
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(line_width);

    auto ac = new VectorActorControler(name + label[dim], actor);
    ac->setColor(colour[dim]);
    insert(ac);
  }

  return true;
}

bool
VtkWrapper::drawGraph(
  std::string name,
  const std::vector<Eigen::Vector3d>& points,
  const std::vector<Eigen::Matrix<IndexType, 2, 1>>& vertices_pairs,
  std::vector<double>* scalars,
  double line_width)
{
  auto data = processPolyData(points, vertices_pairs);

  vtkNew<vtkPolyDataMapper> mapper;

  if (scalars != nullptr) {
    double min = (*scalars)[0];
    double max = min;
    vtkNew<vtkDoubleArray> darray;
    darray->Allocate(scalars->size());
    for (auto d : *scalars) {
      darray->InsertNextValue(d);
      if (d > max)
        max = d;
      if (d < min)
        min = d;
    }

    data->GetCellData()->SetScalars(darray);
    mapper->SetScalarModeToUseCellData();
    mapper->SetColorModeToMapScalars();
    // mapper->SetScalarRange(min, max);
    mapper->ScalarVisibilityOn();
    // add scalarbar
    // add scalarbar
    auto hueLut = get_lookuptable(min, max);
    mapper->SetLookupTable(hueLut);
    vtkNew<vtkScalarBarActor> bar;
    bar->SetLookupTable(mapper->GetLookupTable());
    bar->SetNumberOfLabels(5);
    auto ac2 = new Actor2dControler(name + "_bar", bar.Get());
    insert(ac2);
  } else {
  }

  mapper->SetInputData(data);
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);

  auto ac = new Actor3dControler(name, actor);
  ac->setSize(line_width);
  if (!scalars)
    ac->setIntColor(220, 87, 18);
  insert(ac);

  return true;
}

bool
VtkWrapper::drawGraph(std::string name,
                      const RowMatrix<double>& points,
                      const RowMatrix<IndexType>& edge_pairs,
                      const RowMatrix<double>* scalars,
                      double line_width)
{
  auto data = processPolyData(points, edge_pairs);

  vtkNew<vtkPolyDataMapper> mapper;

  if (scalars != nullptr) {
    size_t n_s = scalars->rows();
    double min = (*scalars)(0, 0);
    double max = min;
    vtkNew<vtkDoubleArray> darray;
    darray->Allocate(scalars->size());
    for (size_t i = 0; i < n_s; i++) {
      double d = (*scalars)(i, 0);
      darray->InsertNextValue(d);
      if (d > max)
        max = d;
      if (d < min)
        min = d;
    }
    data->GetCellData()->SetScalars(darray);
    mapper->SetScalarModeToUseCellData();
    mapper->SetColorModeToMapScalars();
    // mapper->SetColorModeToMapScalars();
    mapper->ScalarVisibilityOn();
    // add scalarbar
    auto hueLut = get_lookuptable(min, max);
    mapper->SetLookupTable(hueLut);
    vtkNew<vtkScalarBarActor> bar;
    bar->SetLookupTable(hueLut);
    // bar->SetLookupTable(mapper->GetLookupTable());
    bar->SetNumberOfLabels(5);
    auto ac2 = new Actor2dControler(name + "_bar", bar.Get());
    insert(ac2);
  } else {
  }

  mapper->SetInputData(data);
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);

  auto ac = new Actor3dControler(name, actor);
  ac->setSize(line_width);
  insert(ac);
  if (!scalars) {
    ac->setIntColor(220, 87, 18);
  }

  return true;
}

bool
VtkWrapper::showIds(std::string name)
{
  auto target = findActor3d(name);
  if (std::get<0>(target)) {
    std::cerr << name + " not found" << std::endl;
    return false;
  }
  auto& ac = std::get<1>(target);
  // only support show mesh's id now
  if (ac->getClassType() != Actor3dControler::ActorType::MESH) {
    std::cerr << name + " is not mesh" << std::endl;
    return false;
  }

  MeshActorControler* mac = static_cast<MeshActorControler*>(ac);
  auto data = mac->get_data();
  vtkNew<vtkIdFilter> id_filter;
  id_filter->SetInputData(data);
  id_filter->PointIdsOn();

  vtkNew<vtkLabeledDataMapper> pLabel;
  pLabel->SetInputConnection(id_filter->GetOutputPort());
  pLabel->SetLabelModeToLabelIds();

  vtkNew<vtkActor2D> pActor;
  pActor->SetMapper(pLabel);

  auto ac2d = new Actor2dControler(name + "_ids", pActor);
  insert(ac2d);

  return true;
}

bool
VtkWrapper::drawIsosurface(std::string name,
                           const RowMatrix<double>& points,
                           const RowMatrix<IndexType>& faces,
                           const RowMatrix<double>& scalars,
                           double value)
{
  auto data = processPolyData(points, faces);
  auto values = convert_arrary(scalars);
  data->GetPointData()->SetScalars(values);
  vtkNew<vtkContourFilter> filter;
  filter->SetInputData(data);
  filter->SetValue(0, value);
  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(filter->GetOutputPort());
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  auto ac = new Actor3dControler(name, actor);
  insert(ac);
  return true;
}

bool
VtkWrapper::drawLabels(std::string name,
                       const RowMatrix<double>& points,
                       const RowMatrix<int>& labels)
{
  assert(points.rows() == labels.rows());

  size_t n = points.rows();

  auto data = processPointData(points);

  vtkNew<vtkIntArray> slabels;
  // labels
  slabels->SetNumberOfValues(n);
  slabels->SetName("labels");
  for (int i = 0; i < n; i++) {
    // slabels->SetValue(i, i);
    slabels->SetValue(i, labels(i));
  }

  auto point_source = data->GetPointData();
  point_source->SetScalars(slabels);

  // render labels

  vtkNew<vtkLabeledDataMapper> mapper;
  mapper->SetInputData(data);
  mapper->SetLabelModeToLabelScalars();
  vtkNew<vtkActor2D> actor2d;
  actor2d->SetMapper(mapper);

  auto ac2d = new Actor2dControler(name, actor2d);
  insert(ac2d);

  return true;
}

bool
VtkWrapper::drawLabels(std::string name,
                       const RowMatrix<double>& points,
                       const std::vector<int>* labels)
{
  size_t n = points.rows();
  // labels
  vtkNew<vtkIntArray> slabels;
  slabels->SetNumberOfValues(n);
  slabels->SetName("labels");

  if (labels) {
    auto& tabs = *labels;
    assert(points.rows() == tabs.size());
    for (int i = 0; i < n; i++) {
      slabels->SetValue(i, tabs[i]);
    }
  } else {
    for (int i = 0; i < n; ++i)
      slabels->SetValue(i, i);
  }

  auto data = processPointData(points);
  auto point_source = data->GetPointData();
  point_source->SetScalars(slabels);

  // render labels

  vtkNew<vtkLabeledDataMapper> mapper;
  mapper->SetInputData(data);
  mapper->SetLabelModeToLabelScalars();
  vtkNew<vtkActor2D> actor2d;
  actor2d->SetMapper(mapper);

  auto ac2d = new Actor2dControler(name, actor2d);
  insert(ac2d);

  return true;
}

bool
VtkWrapper::drawLabelsFloat(std::string name,
                            const RowMatrix<double>& points,
                            const std::vector<double>* labels)
{
  size_t n = points.rows();
  // labels
  vtkNew<vtkDoubleArray> slabels;
  slabels->SetNumberOfValues(n);
  slabels->SetName("labels");

  if (labels) {
    auto& tabs = *labels;
    assert(points.rows() == tabs.size());
    for (int i = 0; i < n; i++) {
      slabels->SetValue(i, tabs[i]);
    }
  } else {
    for (int i = 0; i < n; ++i)
      slabels->SetValue(i, i);
  }

  auto data = processPointData(points);
  auto point_source = data->GetPointData();
  point_source->SetScalars(slabels);

  // render labels

  vtkNew<vtkLabeledDataMapper> mapper;
  mapper->SetInputData(data);
  mapper->SetLabelModeToLabelScalars();
  vtkNew<vtkActor2D> actor2d;
  actor2d->SetMapper(mapper);

  auto ac2d = new Actor2dControler(name, actor2d);
  insert(ac2d);

  return true;
}

bool
VtkWrapper::drawLines(std::string name,
                      const std::vector<std::vector<Eigen::Vector3d>>& points,
                      bool is_loop)
{
  /* insert vertices */
  // number of all points
  size_t n_points = 0;
  // number of all segments
  size_t n_segs = 0;
  for (auto s : points) {
    n_points += s.size();
    n_segs += s.size() - 1;
  }
  vtkNew<vtkPoints> nodes;
  nodes->GetData()->Allocate(n_points);

  int i = 0;
  for (auto s : points) {
    for (auto p : s) {
      nodes->InsertPoint(i, p.data());
      i++;
    }
  }

  /* insert using vtkLine */
  /* insert polys */
  vtkNew<vtkCellArray> cells;

  vtkNew<vtkPolyData> data;
  //每个单元有1+2个数据 1个存储顶点个数cell_n，cell_n个存储面片顶点的标号
  // cells->GetData()->Allocate((1 + 2) * n_segs);
  int p_id = 0;
  int s_id; // if draw loop, we should add the start point at the end
  for (auto s : points) {
    vtkNew<vtkPolyLine> polyline;
    if (is_loop)
      polyline->GetPointIds()->SetNumberOfIds(s.size() + 1);
    else
      polyline->GetPointIds()->SetNumberOfIds(s.size());
    s_id = p_id;
    for (int i = 0; i < s.size(); ++i) {
      // pairs[0] = p_id;
      // pairs[1] = p_id + 1;
      // cells->InsertNextCell(2, pairs);
      polyline->GetPointIds()->SetId(i, p_id);
      p_id++;
    }
    if (is_loop) {
      /*pairs[0] = p_id;
      pairs[1] = p_id - s.size() + 1;
      cells->InsertNextCell(2, pairs);*/
      polyline->GetPointIds()->SetId(s.size(), s_id);
    }
    cells->InsertNextCell(polyline);
    // skip last point of each segment
    // p_id++;
  }
  /* form poly data */

  data->SetPoints(nodes);
  data->SetLines(cells);
  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputData(data);
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->SetLineWidth(2);
  // actor->GetProperty()->SetRepresentationToWireframe();
  // actor->GetProperty()->SetEdgeColor(1.0, 0.4, 0.2);
  auto ac = new MeshActorControler(name, actor);
  // ac->setRenderSyle(1);
  double color[] = { 1.0, 0.4, 0.2 };
  ac->setColor(color);
  insert(ac);
  return true;
}

bool
VtkWrapper::drawMolecule(std::string name,
                         const RowMatrix<double>& points,
                         const RowMatrix<IndexType>& bonds,
                         const RowMatrix<IndexType>* scalars)
{
  size_t n_s = points.rows();

  // auto data = processPolyData(points, bonds);

  // vtkNew<vtkIdTypeArray> darray;
  // darray->Allocate(n_s);
  // if (scalars != nullptr) {
  //  for (size_t i = 0; i < n_s; i++) {
  //    IndexType d = (*scalars)(i, 0);
  //    darray->InsertNextValue(d);
  //  }
  //} else {
  //  for (size_t i = 0; i < n_s; i++) {
  //    darray->InsertNextValue(1);
  //  }
  //}
  // data->GetPointData()->SetScalars(darray);

  //// change polydata to molecule
  // vtkNew<vtkPointSetToMoleculeFilter> filter;
  // filter->SetInputData(data);

  vtkNew<vtkMolecule> mole;
  double const* v;
  for (size_t i = 0; i < n_s; ++i) {
    v = points.row(i).data();
    mole->AppendAtom(1, v[0], v[1], v[2]);
  }

  size_t n_b = bonds.rows();
  for (size_t i = 0; i < n_b; ++i) {
    mole->AppendBond(bonds(i, 0), bonds(i, 1));
  }

  vtkNew<vtkMoleculeMapper> mapper;
  mapper->SetInputData(mole);
  // mapper->SetInputConnection(filter->GetOutputPort());

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);

  auto ac = new Actor3dControler(name, actor);
  insert(ac);

  return true;
}

bool
VtkWrapper::drawMoledule(std::string name,
                         const RowMatrix<double>& points,
                         std::vector<unsigned short>* label)
{
  size_t n_s = points.rows();
  vtkNew<vtkMolecule> mole;
  double const* v;
  auto group_ids = *label;
  if (label) {
    for (size_t i = 0; i < n_s; ++i) {
      v = points.row(i).data();
      mole->AppendAtom(group_ids[i], v[0], v[1], v[2]);
    }
  } else {
    for (size_t i = 0; i < n_s; ++i) {
      v = points.row(i).data();
      mole->AppendAtom(i, v[0], v[1], v[2]);
    }
  }

  // add color
  size_t cn = _preset_colors.size();
  auto mn = mole->GetNumberOfAtoms();
  vtkNew<vtkDoubleArray> colors;
  colors->SetName("Colors");
  colors->SetNumberOfComponents(3);
  colors->Allocate(3 * mn);
  for (IndexType i = 0; i < mn; i++) {
    auto atomic_n = mole->GetAtomAtomicNumber(i);
    colors->InsertNextTypedTuple(_preset_colors[atomic_n % cn].data());
  }
  mole->GetAtomData()->AddArray(colors);

  vtkNew<vtkMoleculeMapper> mapper;
  mapper->SetInputData(mole);
  /*mapper->SetInputArrayToProcess(
      0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_VERTICES, "Colors");*/

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);

  auto ac = new Actor3dControler(name, actor);
  insert(ac);

  return true;
}

bool
VtkWrapper::drawPoints(std::string name,
                       const std::vector<Eigen::Vector3d>& points,
                       double point_size,
                       std::vector<unsigned short>* label)
{
  // std::vector<vtkFacetTuple<1>> polys;
  std::vector<VertexList<1>> polys;
  polys.resize(points.size());
  size_t n = points.size();
  vtkNew<vtkPolyDataMapper> mapper;
  for (int i = 0; i < n; ++i) {
    polys[i] = VertexList<1>(i);
  }

  auto data = processPolyData(points, polys);
  if (label != nullptr) {
    vtkNew<vtkDoubleArray> scalar;
    scalar->Allocate(n);

    vtkNew<vtkLookupTable> lut;
    size_t color_n = _preset_colors.size();
    lut->SetNumberOfTableValues(static_cast<vtkIdType>(color_n));
    lut->Build();
    for (size_t i = 0; i < color_n; i++) {
      lut->SetTableValue(
        i, _preset_colors[i][0], _preset_colors[i][1], _preset_colors[i][2]);
    }

    for (int i = 0; i < n; ++i) {
      scalar->InsertNextValue(static_cast<double>((*label)[i]));
    }
    data->GetPointData()->SetScalars(scalar);
    mapper->SetScalarRange(0, color_n);
    mapper->SetLookupTable(lut);
    mapper->SetScalarModeToUsePointData();
  }
  mapper->SetInputData(data);
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);

  auto ac = new PointsActorControler(name, actor);
  ac->setSize(point_size);
  insert(ac);

  return true;
}

bool
VtkWrapper::drawPoints(std::string name,
                       const RowMatrix<double>& points,
                       double point_size,
                       std::vector<unsigned short>* label)
{
  size_t n = points.rows();
  vtkNew<vtkPolyDataMapper> mapper;
  RowMatrix<IndexType> polys(n, 1);
  for (int i = 0; i < n; ++i)
    polys(i) = i;

  auto data = processPolyData(points, polys);
  if (label != nullptr) {
    vtkNew<vtkDoubleArray> scalar;
    scalar->Allocate(n);

    vtkNew<vtkLookupTable> lut;
    size_t color_n = _preset_colors.size();
    lut->SetNumberOfTableValues(static_cast<vtkIdType>(color_n));
    lut->Build();
    for (size_t i = 0; i < color_n; i++) {
      lut->SetTableValue(
        i, _preset_colors[i][0], _preset_colors[i][1], _preset_colors[i][2]);
    }

    for (int i = 0; i < n; ++i) {
      scalar->InsertNextValue(static_cast<double>((*label)[i]));
    }
    data->GetPointData()->SetScalars(scalar);
    mapper->SetScalarRange(0, color_n);
    mapper->SetLookupTable(lut);
    mapper->SetScalarModeToUsePointData();
  }
  mapper->SetInputData(data);
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);

  auto ac = new PointsActorControler(name, actor);
  ac->setSize(point_size);
  insert(ac);

  return true;
}

bool
VtkWrapper::drawSegments(
  std::string name,
  const std::vector<Eigen::Vector3d>& points,
  const std::vector<Eigen::Matrix<IndexType, 2, 1>>& vertices_pairs,
  double line_width)
{
  auto data = processPolyData(points, vertices_pairs);

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputData(data);
  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  auto ac = new SegmentActorControler(name, actor);
  ac->setColor(0.2, 0.9, 0.0);
  ac->setSize(line_width);
  insert(ac);

  return true;
}

bool
VtkWrapper::drawTetMesh(
  std::string name,
  const std::vector<Eigen::Vector3d>& points,
  const std::vector<Eigen::Matrix<IndexType, 3, 1>>& faces)
{
  drawMesh<3>(name, points, faces);
  return true;
}

bool
VtkWrapper::drawVector(std::string name,
                       const std::vector<Eigen::Vector3d>& points,
                       const std::vector<Eigen::Vector3d>& vectors,
                       double scale_factor,
                       double line_width)
{
  auto locs = convert_points(points);
  auto dirs = convert_vectors(vectors);

  vtkNew<vtkPolyData> ps;
  ps->SetPoints(locs);
  ps->GetPointData()->SetVectors(dirs);

  vtkNew<vtkHedgeHog> hedgehog;
  hedgehog->SetInputData(ps);
  // hedgehog->SetScaleFactor(scale_factor);
  if (std::abs(scale_factor - 1.0) < 1e-7)
    hedgehog->SetVectorModeToUseVector();
  else
    hedgehog->SetScaleFactor(scale_factor);

  // Mapper
  vtkNew<vtkPolyDataMapper> mapper;
  // mapper->SetInputConnection(input);
  mapper->SetInputConnection(hedgehog->GetOutputPort());

  // Actor in scene _actor->SetMapper(mapper)
  vtkNew<vtkActor> actor;

  // Actor in scene
  actor->SetMapper(mapper);
  actor->GetProperty()->SetLineWidth(line_width);

  auto ac = new VectorActorControler(name, actor);
  ac->setColor({ 153 / 255.0, 50 / 255.0, 204 / 255.0 });

  // actor->GetProperty()->SetColor(color.data());

  insert(ac);

  return true;
}

bool
VtkWrapper::drawVector(std::string name,
                       const RowMatrix<double>& points,
                       const RowMatrix<double>& vectors,
                       double scale_factor,
                       double line_width)
{
  auto locs = convert_points(points);
  auto dirs = convert_vectors(vectors);

  vtkNew<vtkPolyData> pd;
  pd->SetPoints(locs);
  pd->GetPointData()->SetVectors(dirs);

  vtkNew<vtkHedgeHog> hedgehog;
  hedgehog->SetInputData(pd);
  if (std::abs(scale_factor - 1.0) < 1e-7)
    hedgehog->SetVectorModeToUseVector();
  else
    hedgehog->SetScaleFactor(scale_factor);

  // Mapper
  vtkNew<vtkPolyDataMapper> mapper;
  // mapper->SetInputConnection(input);
  mapper->SetInputConnection(hedgehog->GetOutputPort());

  // Actor in scene _actor->SetMapper(mapper)
  vtkNew<vtkActor> actor;

  // Actor in scene
  actor->SetMapper(mapper);
  actor->GetProperty()->SetLineWidth(line_width);

  auto ac = new VectorActorControler(name, actor);
  ac->setColor({ 153 / 255.0, 50 / 255.0, 204 / 255.0 });

  // actor->GetProperty()->SetColor(color.data());

  insert(ac);
  return true;
}

bool
VtkWrapper::drawVectorWithScalars(std::string name,
                                  const std::vector<Eigen::Vector3d>& points,
                                  const std::vector<Eigen::Vector3d>& vectors,
                                  std::vector<double>& scalars,
                                  double scale_factor,
                                  double line_width)
{
  size_t n = points.size();
  auto locs = convert_points(points);
  auto dirs = convert_vectors(vectors);
  auto scals = convert_arrary(scalars);

  double low_bound = 1e20, high_bound = -1e20;
  for (int i = 0; i < n; ++i) {
    if (scalars[i] < low_bound) {
      low_bound = scalars[i];
    } else if (scalars[i] > high_bound) {
      high_bound = scalars[i];
    }
  }

  // std::cout << "stress low : " << low_bound << "--> high : " << high_bound
  //           << std::endl;

  vtkNew<vtkPolyData> ps;
  ps->SetPoints(locs);
  ps->GetPointData()->SetVectors(dirs);
  ps->GetPointData()->SetScalars(scals);

  vtkNew<vtkHedgeHog> hedgehog;
  hedgehog->SetInputData(ps);
  if (std::abs(scale_factor - 1.0) < 1e-7)
    hedgehog->SetVectorModeToUseVector();
  else
    hedgehog->SetScaleFactor(scale_factor);

  // Mapper
  vtkNew<vtkPolyDataMapper> mapper;
  // mapper->SetInputConnection(input);
  mapper->SetInputConnection(hedgehog->GetOutputPort());
  mapper->SetScalarModeToUsePointData();
  mapper->SetScalarRange(low_bound, high_bound);

  // Actor in scene _actor->SetMapper(mapper)
  vtkNew<vtkActor> actor;

  // Actor in scene
  actor->SetMapper(mapper);
  actor->GetProperty()->SetLineWidth(line_width);

  auto ac = new VectorActorControler(name, actor);

  insert(ac);

  // actor->GetProperty()->SetColor(color.data());

  return true;
}

void
VtkWrapper::setOpacity(std::string name, double opacity)
{
  auto target = findActor3d(name);
  if (std::get<0>(target)) {
    std::get<1>(target)->setOpacity(opacity);
  }
}

void
VtkWrapper::setVisibility(std::string name, bool flag)
{
  auto target = findActor3d(name);
  if (std::get<0>(target)) {
    std::get<1>(target)->setVisibility(flag);
  } else {
    auto t2 = findActor2d(name);
    if (std::get<0>(t2)) {
      std::get<1>(t2)->setVisibility(flag);
    }
  }
}

void
VtkWrapper::setColor(std::string name, double* color)
{
  auto target = findActor3d(name);
  if (std::get<0>(target)) {
    std::get<1>(target)->setColor(color);
  }
}

std::tuple<int, int, int>
VtkWrapper::getColor(std::string name)
{
  auto target = findActor3d(name);
  if (std::get<0>(target)) {
    return std::get<1>(target)->getIntColor();
  }
  return std::tuple<int, int, int>(0, 0, 0);
}

void
VtkWrapper::setSize(std::string name, double size)
{
  auto target = findActor3d(name);
  if (std::get<0>(target)) {
    std::get<1>(target)->setSize(size);
  }
}

double
VtkWrapper::getSize(std::string name)
{
  auto target = findActor3d(name);
  if (std::get<0>(target)) {
    return std::get<1>(target)->getSize();
  }
  return 0.0;
}

std::tuple<double, double>
VtkWrapper::getStatus(std::string name)
{
  auto target = findActor3d(name);
  if (std::get<0>(target)) {
    auto status = std::get<1>(target)->getStatus();
    return std::tuple<double, double>(std::get<1>(status), std::get<2>(status));
  }
  return std::tuple<double, double>(1.0, 0.0);
}

void
VtkWrapper::setRenderStyle(std::string name, int render_style)
{
  auto target = findActor3d(name);
  if (std::get<0>(target)) {
    std::get<1>(target)->setRenderSyle(render_style);
  }
}

void
VtkWrapper::setVertexScalars(std::string name,
                             std::vector<double>& scalars,
                             double lower_bound,
                             double upper_bound)
{
  auto target = findActor3d(name);
  if (std::get<0>(target))
    return;
  auto& ac = std::get<1>(target);
  auto actor = ac->get_actor();
  auto mapper = actor->GetMapper();
  auto data = mapper->GetInput();
  vtkNew<vtkDoubleArray> darray;
  darray->Allocate(scalars.size());
  for (auto d : scalars) {
    darray->InsertNextValue(d);
  }
  darray->SetName(name.c_str());
  data->GetPointData()->SetScalars(darray);
  // data->GetPointData()->SetActiveScalars("scalars");

  mapper->SetScalarModeToUsePointData();
  mapper->SetScalarRange(lower_bound, upper_bound);

  refresh();
}

bool
VtkWrapper::exist(std::string name)
{
  auto target = findActor3d(name);
  return std::get<0>(target);
}

void
VtkWrapper::testRenderFunction()
{}

} // namespace viewtools