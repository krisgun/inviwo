/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2013-2015 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 *********************************************************************************/

#ifndef IVW_CANVASQT_H
#define IVW_CANVASQT_H

#include <modules/openglqt/openglqtmoduledefine.h>
#include <modules/openglqt/canvaseventmapper.h>
#include <modules/opengl/canvasgl.h>
#include <inviwo/qt/widgets/eventconverterqt.h>
#include <inviwo/core/network/processornetworkevaluator.h>
#include <inviwo/core/common/inviwo.h>

#define QT_NO_OPENGL_ES_2
#define GLEXT_64_TYPES_DEFINED

#include <warn/push>
#include <warn/ignore/all>
#include <QGLWidget>
#include <warn/pop>

namespace inviwo {

class IVW_MODULE_OPENGLQT_API CanvasQt : public QGLWidget, public CanvasGL {
    friend class CanvasProcessorWidgetQt;

public:
    explicit CanvasQt(QGLWidget* parent = nullptr, uvec2 dim = uvec2(256,256));
    ~CanvasQt() = default;

    static void defineDefaultContextFormat();

    virtual void activate() override;
    virtual void glSwapBuffers() override;
    virtual void update() override;
    void repaint();

    virtual void resize(uvec2 size) override;

    virtual std::unique_ptr<Canvas> create() override;

protected:
    void initializeGL() override;
    void paintGL() override;

    virtual bool event(QEvent *e) override;

    virtual void resizeEvent(QResizeEvent* event) override;
    static CanvasQt* getSharedCanvas();

private:
    void touchFallback(QTouchEvent*);


    static QGLWidget* sharedGLContext_; //For rendering-context sharing

    static CanvasQt* sharedCanvas_;
    static QGLFormat sharedFormat_;
    bool swapBuffersAllowed_;
    CanvasEventMapper eventMapper_;
};

} // namespace

#endif // IVW_CANVASQT_H