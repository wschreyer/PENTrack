// Copyright (c) 2018  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.1/Three/include/CGAL/Three/Scene_item_rendering_helper.h $
// $Id: Scene_item_rendering_helper.h 9072324 %aI Maxime Gimeno
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Maxime Gimeno

#ifndef SCENE_ITEM_RENDERING_HELPER_H
#define SCENE_ITEM_RENDERING_HELPER_H

#include <CGAL/license/Three.h>


#include <CGAL/Three/Scene_item.h>
#include <QRunnable>


#ifdef demo_framework_EXPORTS
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_EXPORT
#else
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_IMPORT
#endif

struct PRIV;
class QMenu;
class QSlider;
namespace CGAL {
namespace Three{

class Viewer_interface;
struct Triangle_container;
struct Edge_container;
struct Point_container;

//!
//! \brief The Scene_item_rendering_helper class is a convenience class for  constructing an item.
//! It is more elaborated than a `Scene_item` and facilitates the process of creating an item that can
//! be rendered.
//!
class DEMO_FRAMEWORK_EXPORT Scene_item_rendering_helper
    :public Scene_item
{
  Q_OBJECT
public:
  Scene_item_rendering_helper();
  ~Scene_item_rendering_helper();

  //!
  //! \brief The `Gl_data_name` enum is used as a flag to specify what should be
  //! re-computed during `computeElements()`. The flag corresponding to this enum is
  //! `Gl_data_names`, and multiple flags can be combined whith the operator `|`.
  //! For instance, you can use `GEOMETRY|COLORS` as a single value.
  //! @todo Review Laurent Rineau We need to find a better name. 1. Do not refer to OpenGL. 2. Why "name"?
  //!
  enum Gl_data_name{
    GEOMETRY = 0x1,                     //!< Invalidates the vertices, edges and faces.
    COLORS   = 0x2,                     //!< Invalidates the color of each vertex
    NORMALS  = 0x4,                     //!< Invalidate the normal of each vertex.
    ALL      = GEOMETRY|COLORS|NORMALS  //!< Invalidate everything
  };
#ifdef DOXYGEN_RUNNING
  //! \brief Flag interface for Scene_item::Gl_data_name.
  //! \todo Review Laurent Rineau:  Should be explained better. Points to `QFlags`...
  enum Gl_data_names{};
#endif
  Q_DECLARE_FLAGS(Gl_data_names, Gl_data_name)

  QMenu* contextMenu() Q_DECL_OVERRIDE;

  /*!
     * \brief processData calls `computeElements()` 
     * 
     * @todo in a dedicated thread so the
     * application does not get stuck while the processing is performed.
     * Emits `dataProcessed()`.
     */
   virtual void processData(Gl_data_names name) const;
   
  //!
  //! \brief setAlpha sets the integer value of the alpha channel of this item.
  //! Also updates the slider value.
  //! It must be between 0 and 255.
  //! \param alpha the integer value for the alpha channel.
  //!
  void setAlpha(int alpha) Q_DECL_OVERRIDE;
  //! \brief The item's bounding box.
  //!
  //! If the Bbox has never been computed, computes it and
  //! saves the result for further calls.
  //! @returns the item's bounding box.
  Scene_item::Bbox bbox()const Q_DECL_OVERRIDE;
  //!
  //! \brief getTriangleContainer returns the `id`th `Triangle_container`.
  //!
  CGAL::Three::Triangle_container* getTriangleContainer(std::size_t id) const;
  //!
  //! \brief getEdgeContainer returns the `id`th `Edge_container`.
  //!
  CGAL::Three::Edge_container* getEdgeContainer(std::size_t id)const;
  //!
  //! \brief getPointContainer returns the `id`th `Point_container`.
  //!
  CGAL::Three::Point_container* getPointContainer(std::size_t id)const;

  //!
  //! \brief setTriangleContainer sets the `id`th `Triangle_container` to `tc`.
  //!
  //! If `id` is bigger than the current size of the container vector, this vector is
  //! resized accordingly. This means that for optimisation reasons, containers should be created
  //! decreasingly.
  //!
  void setTriangleContainer(std::size_t id,
                            Triangle_container* tc);

  //!
  //! \brief setEdgeContainer sets the `id`th `Edge_container` to `tc`.
  //!
  //! If `id` is bigger than the current size of the container vector, this vector is
  //! resized accordingly. This means that for optimisation reasons, containers should be created
  //! decreasingly.
  //!
  void setEdgeContainer(std::size_t id,
                        Edge_container* tc);
  
  //!
  //! \brief setPointContainer sets the `id`th `Point_container` to `tc`.
  //!
  //! If `id` is bigger than the current size of the container vector, this vector is
  //! resized accordingly. This means that for optimisation reasons, containers should be created
  //! decreasingly.
  //!
  void setPointContainer(std::size_t id,
                        Point_container* tc);

  //!
  //! \brief setBuffersFilled specifies if the data should be re-computed.
  //!
  //! If called with `false`, the item rendering data will be re-computed at the next `draw()`.
  //! If called with `true`, the item rendering data is considered ready and will not be computed
  //! until `setBuffersFilled()` is called with `false` again.
  //!
  void setBuffersFilled(bool b) const;

  //!
  //! \brief getBuffersFilled returns `false` if the item rendering data needs to be re-computed.,
  //! `true` otherwise.
  //! \see `setBuffersFilled()`
  bool getBuffersFilled()const;

  //!
  //! \brief getBuffersInit returns true if the `Vao`s of `viewer` are ready
  //! for rendering.
  //!
  bool getBuffersInit(Viewer_interface *viewer)const;

  //!
  //! \brief setBuffersInit specifies if the `Vbo`s need to be initialized.
  //!
  //! If called with `false`, the item `Vbo`s will be refilled at the next `draw()`.
  //! If called with `true`, the item `Vbo`s are considered ready and will not be refilled
  //! until `setBuffersInit()` is called with `false` again.
  //!
  //! This function should be called in the drawing functions, when `getBuffersFilled()` is `true`.
  //!
  void setBuffersInit(Viewer_interface *viewer, bool val) const;
  
  //! \brief the item's bounding box's diagonal length.
  //!
  //! If the diagonal's length has never been computed, computes it and
  //! saves the result for further calls.
  //! @returns the item's bounding box's diagonal length.
  //! @todo must replace the one from Scene_item eventually
  virtual double diagonalBbox() const Q_DECL_OVERRIDE;
protected:


  //!Returns a pointer to the slider initialized in initGL();
  QSlider* alphaSlider();

  //!Returns`true` if `initGL()` was called.
  bool isInit()const;

  //!Returns the float alpha value of an item.
  //! This value is between 0.0f and 1.0f.
  float alpha() const Q_DECL_OVERRIDE;

  /*! Fills the `Vbo`s with data.
     */
  virtual void initializeBuffers(Viewer_interface*)const{}

  //!Creates the VAOs and VBOs for each existing viewer.
  virtual void initGL() const;
  //!
  //! Computes the items Bbox and stores the result. Must be overridden.
  //!@todo must replace the one from Scene_item eventually.
  virtual void compute_bbox() const Q_DECL_OVERRIDE = 0;
  //!
  //! \brief setBbox allows to set the Bbox in compute_bbox();
  //! \param b
  //!
  void setBbox(Bbox b);
  
  virtual void computeElements()const{}
private:
  friend struct PRIV;
  mutable PRIV* priv;
};//end Scene_item_rendering_helper

}}

#endif // SCENE_ITEM_RENDERING_HELPER_H
