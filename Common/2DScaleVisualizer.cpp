#include "stdafx.h"
#include "2DScaleVisualizer.h"

//------------------- class GLViewport ------------------

void GLViewport::Set(GLint x, GLint y, GLsizei width, GLsizei height) {
	m_X = x;
	m_Y = y;
	m_Width = width;
	m_Height = height;
}

void GLViewport::Set(const GLint *viewport) {
	m_X = viewport[0];
	m_Y = viewport[1];
	m_Width = viewport[2];
	m_Height = viewport[3];
}

void GLViewport::Set(const CRect &rect) {
	m_X = rect.left;
	m_Y = rect.bottom;
	m_Width = rect.Width();
	m_Height = rect.Height();
}

void GLViewport::ShiftBy(GLint x, GLint y) {
	m_X += x;
	m_Y += y;
}

void GLViewport::ShiftBy(CPoint point) {
	m_X += point.x;
	m_Y += point.y;
}

void GLViewport::ResizeLeftBottom(GLint x, GLint y) {
	m_X += x;
	m_Y += y;
	m_Width -= x;
	m_Height += y;
}

void GLViewport::ResizeRightTop(GLint width, GLint height) {
	m_Width += width;
	m_Height -= height;
}

void GLViewport::FillArray(GLint out[4]) const {
	out[0] = m_X;
	out[1] = m_Y;
	out[2] = m_Width;
	out[3] = m_Height;
}

/**/
//------------------- class ScaleVisualizer --------------
/**/

CScaleVisualizer::CScaleVisualizer()
{
	m_IsMoved = false;
	m_Trasparency = 0.6;
	m_DisplayListIndex = 0;
	m_ScaleTop = 0.9;
	m_ScaleBottom = -0.9;
	m_ScaleLeft = 0.7;
	m_ScaleRight = 1;
	m_FontHeight = 0.05;

	m_Font = nullptr;
	//m_ClientGLintArray = {0, 0, 0, 0}; //TODO: инициализацию бы
}

void CScaleVisualizer::SetTrasparency(int percent) {
	m_Trasparency = percent / 100.0;

	ClearDisplayList();
}

void CScaleVisualizer::DrawScale(const C2DScale &scale){
	
		SetProjection();

		if (m_DisplayListIndex == 0)
		{
			bool isDescrete = scale.IsDescrete();
			size_t count = scale.GetCount();
			if (count == 1) {
				count++;
			}

			//не делим на 0
			double factor = (m_ScaleTop - m_ScaleBottom) / (count - 1);

			m_DisplayListIndex = glGenLists(1);
			glNewList(m_DisplayListIndex, GL_COMPILE);
			for (size_t i = 0; i < count - 1; i++)
			{
				// Квадратик с цветом шкалы
				C2DColor color = scale.GetColor(i);
				glColor4d(color.GetR(), color.GetG(), color.GetB(), m_Trasparency);
				glBegin(GL_QUADS);
				glVertex2d(m_ScaleLeft, m_ScaleBottom + factor * i);
				glVertex2d(m_ScaleRight, m_ScaleBottom + factor * i);

				if (!isDescrete) {
					C2DColor topColor = scale.GetColor(i + 1);
					glColor4d(topColor.GetR(), topColor.GetG(), topColor.GetB(), m_Trasparency);
				}

				glVertex2d(m_ScaleRight, m_ScaleBottom + factor * (i + 1));
				glVertex2d(m_ScaleLeft, m_ScaleBottom + factor * (i + 1));

				glEnd();
			}

			/*for (UNINT i = 0; i < count; i++){
			DrawScaleMark(m_ScaleBottom + factor * i, scale.GetStringValue(i), isDescrete);*/
			//}

			//// Чёрная линия над верхним квадратиком
			//glColor4d(0, 0, 0, m_Trasparency);
			//	glBegin(GL_LINES);
			//		glVertex2d(m_ScaleLeft, m_ScaleTop);
			//		glVertex2d(m_ScaleRight, m_ScaleTop);
			//	glEnd();				

			DrawScaleMarks(scale);

			// Рамочка вокруг шкалы
			DrawScaleBorder();

			// Рамочка вокруг Viewport
			DrawViewportBorder();

			glEndList();
		}

		//glColor3d(0, 0, 0);
		glCallList(m_DisplayListIndex);

		RestoreProjection();
	}

void CScaleVisualizer::ClearDisplayList() {
	if (m_DisplayListIndex != 0){
		glDeleteLists(m_DisplayListIndex, 1);
		m_DisplayListIndex = 0;
	}
}

void CScaleVisualizer::SetProjection(){
	glGetIntegerv(GL_VIEWPORT, m_ClientGLintArray);
	GLint scaleGLintArray[4];
	
	m_ScaleViewport.FillArray(scaleGLintArray);
	m_FontHeight = 25.0 / scaleGLintArray[3];
	glViewport(scaleGLintArray[0], m_ClientGLintArray[3] - scaleGLintArray[1], scaleGLintArray[2], scaleGLintArray[3]);
	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
}

void CScaleVisualizer::RestoreProjection(){
	glViewport(m_ClientGLintArray[0], m_ClientGLintArray[1], m_ClientGLintArray[2], m_ClientGLintArray[3]);
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

void CScaleVisualizer::DrawViewportBorder()
{
	glColor4d(0, 0, 0, m_Trasparency);

	// Линия справа
	glBegin(GL_LINES);
	glVertex2d(1, -1);
	glVertex2d(1, 1);
	glEnd();

	// Линия сверху
	glBegin(GL_LINES);
	glVertex2d(1, 1);
	glVertex2d(-1, 1);
	glEnd();

	// Линия слева
	glBegin(GL_LINES);
	glVertex2d(-1, 1);
	glVertex2d(-1, -1);
	glEnd();

	// Линия снизу
	glBegin(GL_LINES);
	glVertex2d(-1, -1);
	glVertex2d(1, -1);
	glEnd();
}

void CScaleVisualizer::DrawScaleBorder(){
	glColor4d(0, 0, 0, m_Trasparency);

	// Линия справа
	glBegin(GL_LINES);
	glVertex2d(m_ScaleRight, m_ScaleBottom);
	glVertex2d(m_ScaleRight, m_ScaleTop);
	glEnd();

	// Линия сверху
	glBegin(GL_LINES);
	glVertex2d(m_ScaleRight, m_ScaleTop);
	glVertex2d(m_ScaleLeft, m_ScaleTop);
	glEnd();

	// Линия слева
	glBegin(GL_LINES);
	glVertex2d(m_ScaleLeft, m_ScaleTop);
	glVertex2d(m_ScaleLeft, m_ScaleBottom);
	glEnd();

	// Линия снизу
	glBegin(GL_LINES);
	glVertex2d(m_ScaleLeft, m_ScaleBottom);
	glVertex2d(m_ScaleRight, m_ScaleBottom);
	glEnd();
}

std::vector<size_t> CScaleVisualizer::GetDrawableMarks(const C2DScale &scale)
{
	std::vector<size_t> result;
	size_t count = scale.GetCount();
	if (count == 1) {
		count++;
	}

	//не делим на 0
	double factor = (m_ScaleTop - m_ScaleBottom) / (count - 1);

	// Самая нижняя отметка
	result.push_back(0);

	for (size_t i = 1; i < count - 1; i++) {
		double y = (m_ScaleBottom + factor * i);

		// Если отметка не налезает на предыдущую, то добавим её
		if (y >(m_ScaleBottom + factor * result.back()) + m_FontHeight) {
			result.push_back(i);
		}
	}

	//if (!scale.IsDescrete()){
	// Удалим последний добавленный элемент, если он залезает на самую верхнюю отметку
	if ((m_ScaleBottom + factor * result.back()) + m_FontHeight > (m_ScaleBottom + factor * (count - 1))) {
		result.pop_back();
	}

	// Самая верхняя отметка
	result.push_back(count - 1);
	//}

	return result;
}

void CScaleVisualizer::DrawScaleMarks(const C2DScale &scale){
	
	glColor4d(0, 0, 0, m_Trasparency);

	//bool isDescrete = scale.IsDescrete();
	size_t count = scale.GetCount();
	if (count == 1) {
		count++;
	}

	//// В дискретной шкале цветов меньше, чем отметок, т.к. самое большое значение
	//// входит в последний цвет
	//if (scale.IsDescrete())
	//	count++;

	//не делим на 0
	double factor = (m_ScaleTop - m_ScaleBottom) / (count - 1);

	for (size_t i = 0; i < count - 1; i++){
		glBegin(GL_LINES);
		glVertex2d(m_ScaleLeft, m_ScaleBottom + factor * i);
		glVertex2d(m_ScaleRight, m_ScaleBottom + factor * i);
		glEnd();

		glBegin(GL_LINES);
		glVertex2d(m_ScaleRight, m_ScaleBottom + factor * (i + 1));
		glVertex2d(m_ScaleLeft, m_ScaleBottom + factor * (i + 1));
		glEnd();
	}

	std::vector<size_t> marksToDraw = GetDrawableMarks(scale);

	//if (isDescrete){
	for (size_t i = 0; i < marksToDraw.size(); i++){
		double y = 0;
		size_t it = marksToDraw[i];

		//if (isDescrete){
		//	y = ((m_ScaleBottom + factor * it) + (m_ScaleBottom + factor * (it + 1))) / 2;
		//}
		//else{
		y = (m_ScaleBottom + factor * it);
		//}


		DrawScaleMarkValue(y, scale.GetStringValue(it, true));
	}
}


void CScaleVisualizer::DrawScaleMarkValue(double y, const CString &value){
	glColor4d(0, 0, 0, m_Trasparency);

	//glBegin(GL_LINES);
	//	glVertex2d(-0.5, y);
	//	glVertex2d(m_ScaleLeft, y);
	//glEnd();


	wchar_t wChar[255];
	wcscpy_s(wChar, value);

	double xStart = m_ScaleLeft - m_Font->GetTextWidth(wChar, m_FontHeight);
	// Чтобы центр высоты текста и линия шкалы совпадали
	double yStart = y - (m_FontHeight / 2);

	//m_Font->OutText2D(wChar, -0.9, y, m_FontHeight);
	m_Font->OutText2D(wChar, xStart, yStart, m_FontHeight);
}