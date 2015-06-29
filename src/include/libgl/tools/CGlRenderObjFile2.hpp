/*
 * Copyright 2010 Martin Schreiber
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef CGL_RENDER_OBJFILE2_HPP
#define CGL_RENDER_OBJFILE2_HPP

#include "lib/CObjFile.hpp"
#include "libgl/CGlMaterial.hpp"
#include "shaders/shader_blinn/CShaderBlinn.hpp"
#include "shaders/shader_blinn_shadow_map/CShaderBlinnShadowMap.hpp"
#include "shaders/shader_blinn_shadow_and_caustic_map/CShaderBlinnShadowAndCausticMap.hpp"

/**
 * \brief render a object file loaded by CObjFile.hpp
 */
class CGlRenderObjFile2
{
	CShaderBlinn shaderBlinn;
	CShaderBlinnShadowMap shaderBlinnShadowMap;
	CShaderBlinnShadowAndCausticMap shaderBlinnShadowAndCausticMap;

	class CMesh
	{
	public:
		CGlMaterial *material;	///< pointer to existing material
		GLuint vertex_start_id;	///< index of first vertex to render
		GLuint vertices_count;	///< number of vertices to render
	};

	class CGroup
	{
	public:
		std::string name;

		GLuint vertex_start_id;
		GLuint vertices_count;	///< number of vertices in group

		int meshes_count;
		CMesh *meshes;

		CGroup()	:
			vertex_start_id(0),
			vertices_count(0),
			meshes_count(0)
		{
			meshes = nullptr;
		}

		~CGroup()
		{
			if (meshes)
				delete meshes;
		}
	};


	int materials_counter;
	CGlMaterial *materials;

	int groups_count;
	CGroup *groups;

	GLfloat *buf_vertices_start_ptr;
	GLfloat *buf_normals_start_ptr;
	GLfloat *buf_texture_coords_start_ptr;

	size_t buffer_size;
	CGlBuffer buffer;

public:
	CError error;				///< error handler
	bool texture_valid;

public:
	/**
	 * load and initialize buffers with object file handler
	 */
	void load(	CObjFile &cObjFile)
	{
		cleanup();

		/*
		 * LOAD MATERIALS
		 */
		materials = new CGlMaterial[cObjFile.materials.size()];

		CGlMaterial *material = materials;
		for(std::list<CObjFile::CMaterial>::iterator i = cObjFile.materials.begin(); i != cObjFile.materials.end(); i++)
		{
			CObjFile::CMaterial &obj_material = *i;

			material->name = obj_material.name;
			material->ambient_color3 = obj_material.ambient_color3;
			material->diffuse_color3 = obj_material.diffuse_color3;
			material->specular_color3 = obj_material.specular_color3;
			material->specular_exponent = obj_material.specular_exponent;

			if (!obj_material.texture_file.empty())
			{
				material->texture0 = new CGlTexture;
				material->texture0->loadFromFile((std::string("data/textures/obj_textures/")+obj_material.texture_file).c_str());
				CError_PtrAppendReturn(material->texture0);
				material->texture0->bind();
				material->texture0->setParam(GL_TEXTURE_WRAP_R, GL_REPEAT);
				material->texture0->setParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
				material->texture0->setParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
				material->texture0->unbind();
			}

			material++;
		}

		/*
		 * LOAD GROUPS AND MESHES
		 */
		int vertices_buffer_size = 0;
		int normals_buffer_size = 0;
		int texture_coords_buffer_size = 0;
		for(std::list<CObjFile::CGroup>::iterator i = cObjFile.groups.begin(); i != cObjFile.groups.end(); i++)
		{
			CObjFile::CGroup &g = *i;
			vertices_buffer_size += g.faces3_count*3*3*sizeof(GLfloat);
			normals_buffer_size += g.faces3_count*3*3*sizeof(GLfloat);
			// TODO: dont allocate texture buffer space if no texture coordinates are used in the whole .obj file
			texture_coords_buffer_size += g.faces3_count*2*3*sizeof(GLfloat);	// !!!allocate texture buffer - independent of it's usage!!!
		}

		buffer_size = vertices_buffer_size + normals_buffer_size + texture_coords_buffer_size;


		buffer.bind();
		buffer.data(buffer_size, nullptr);

		groups_count = cObjFile.groups.size();
		groups = new CGroup[groups_count];

		GLintptr buffer_vertex_offset = 0;
		GLintptr buffer_normal_offset = vertices_buffer_size;
		GLintptr buffer_texture_coord_offset = vertices_buffer_size+normals_buffer_size;

		buf_vertices_start_ptr = (GLfloat*)buffer_vertex_offset;
		buf_normals_start_ptr = (GLfloat*)buffer_normal_offset;
		buf_texture_coords_start_ptr = (GLfloat*)buffer_texture_coord_offset;

		CGroup *group = groups;
		GLuint vertex_start_id = 0;
		for(std::list<CObjFile::CGroup>::iterator gi = cObjFile.groups.begin(); gi != cObjFile.groups.end(); gi++)
		{
			CObjFile::CGroup &obj_group = *gi;

			group->name = obj_group.name;

			/**
			 * load whole vertex/normal/texture data for current group
			 */

			buffer.subData(buffer_vertex_offset, obj_group.faces3_count*3*3*sizeof(GLfloat), obj_group.lin_vertices);
			buffer_vertex_offset += obj_group.faces3_count*3*3*sizeof(GLfloat);

			buffer.subData(buffer_normal_offset, obj_group.faces3_count*3*3*sizeof(GLfloat), obj_group.lin_normals);
			buffer_normal_offset += obj_group.faces3_count*3*3*sizeof(GLfloat);

			if (obj_group.lin_texture_coords != nullptr)
				buffer.subData(buffer_texture_coord_offset, obj_group.faces3_count*2*3*sizeof(GLfloat), obj_group.lin_texture_coords);
			buffer_texture_coord_offset += obj_group.faces3_count*2*3*sizeof(GLfloat);

			CGlErrorCheck();

			group->vertices_count = obj_group.lin_faces3*3;

			group->meshes_count = obj_group.meshes.size();
			group->meshes = new CMesh[group->meshes_count];

			CMesh *mesh = group->meshes;

			for(std::list<CObjFile::CGroupMaterialMesh>::iterator mi = obj_group.meshes.begin(); mi != obj_group.meshes.end(); mi++)
			{
				CObjFile::CGroupMaterialMesh &obj_mesh = *mi;

				mesh->vertex_start_id = vertex_start_id;
				mesh->vertices_count = obj_mesh.vertex_count;

				// search for according material
				int counter = 0;
				std::list<CObjFile::CMaterial>::iterator mmi;
				for(mmi = cObjFile.materials.begin(); mmi != cObjFile.materials.end(); mmi++)
				{
					if (&*mmi == obj_mesh.material)
						break;
					counter++;
				}

				if (mmi == cObjFile.materials.end())
				{
					std::cerr << CError_Location << "material not found: " << obj_mesh.material->name << std::endl;
					return;
				}

				mesh->material = &materials[counter];

				vertex_start_id += mesh->vertices_count;
				mesh++;
			}

			group++;
		}

		buffer.unbind();
	}

	void init()
	{
		groups = nullptr;
		materials = nullptr;

		CError_AppendReturn(shaderBlinn);
		CError_AppendReturn(shaderBlinnShadowMap);
		CError_AppendReturn(shaderBlinnShadowAndCausticMap);
	}

	CGlRenderObjFile2()	:
		buffer(GL_ARRAY_BUFFER)
	{
		init();
		cleanup();
	}


	/**
	 * load and initialize buffers to render an already loaded .obj file
	 */
	CGlRenderObjFile2(CObjFile &cObjFile)	:
			buffer(GL_ARRAY_BUFFER)
	{
		init();
		cleanup();
		load(cObjFile);
	}


	/**
	 * render the meshes from the object file with materials
	 */
	void render(	CGlLights &cGlLights,
					const GLSL::mat4 &pvm_matrix,
					const GLSL::mat4 &view_model_matrix,
					const GLSL::mat3 &view_model_normal_matrix3
				)
	{
		shaderBlinn.use();

		shaderBlinn.setupUniformsMatrices(pvm_matrix, view_model_matrix, view_model_normal_matrix3);
		shaderBlinn.setupUniformsLights(cGlLights, view_model_matrix*cGlLights.light0_world_pos4);

		buffer.bind();

		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, buf_vertices_start_ptr);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, buf_normals_start_ptr);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, buf_texture_coords_start_ptr);
		glEnableVertexAttribArray(2);

		for (int gi = 0; gi < groups_count; gi++)
		{
			CGroup &g = groups[gi];

			for (int mi = 0; mi < g.meshes_count; mi++)
			{
				CMesh &m = g.meshes[mi];

				shaderBlinn.setupUniformsMaterial(*(m.material));

				if (m.material->texture0 != nullptr)
					m.material->texture0->bind();

				glDrawArrays(GL_TRIANGLES, m.vertex_start_id, m.vertices_count);

				if (m.material->texture0 != nullptr)
					m.material->texture0->unbind();
			}
		}


		glDisableVertexAttribArray(2);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(0);

		buffer.unbind();

		shaderBlinn.disable();

		CGlErrorCheck();
	}



	/**
	 * render the meshes from the object file with materials
	 */
	void renderWithoutProgram()
	{
		buffer.bind();

		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, buf_vertices_start_ptr);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, buf_normals_start_ptr);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, buf_texture_coords_start_ptr);
		glEnableVertexAttribArray(2);

		for (int gi = 0; gi < groups_count; gi++)
		{
			CGroup &g = groups[gi];

			for (int mi = 0; mi < g.meshes_count; mi++)
			{
				CMesh &m = g.meshes[mi];

				if (m.material->texture0 != nullptr)
					m.material->texture0->bind();

				glDrawArrays(GL_TRIANGLES, m.vertex_start_id, m.vertices_count);

				if (m.material->texture0 != nullptr)
					m.material->texture0->unbind();
			}
		}


		glDisableVertexAttribArray(2);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(0);

		buffer.unbind();

		CGlErrorCheck();
	}


	/**
	 * render the meshes from the object file using a shadow map texture
	 */
	void renderWithShadowMap(	CGlLights &cGlLights,
								const GLSL::mat4 &pvm_matrix,
								const GLSL::mat4 &view_model_matrix,
								const GLSL::mat3 &view_model_normal_matrix3,
								const GLSL::mat4 &view_matrix,
								CGlTexture &depth_texture,
								const GLSL::mat4 &shadow_map_matrix
				)
	{
		shaderBlinnShadowMap.use();

		shaderBlinnShadowMap.setupUniformsMatrices(pvm_matrix, view_model_matrix, view_model_normal_matrix3);
		shaderBlinnShadowMap.setupUniformsShadowMapping(shadow_map_matrix);
		shaderBlinnShadowMap.setupUniformsLights(cGlLights, view_matrix*cGlLights.light0_world_pos4);

		glActiveTexture(GL_TEXTURE1);
		depth_texture.bind();
		depth_texture.setParam(GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
		depth_texture.setParam(GL_TEXTURE_COMPARE_FUNC, GL_GREATER);

		glActiveTexture(GL_TEXTURE0);

		buffer.bind();

		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, buf_vertices_start_ptr);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, buf_normals_start_ptr);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, buf_texture_coords_start_ptr);
		glEnableVertexAttribArray(2);

		for (int gi = 0; gi < groups_count; gi++)
		{
			CGroup &g = groups[gi];

			for (int mi = 0; mi < g.meshes_count; mi++)
			{
				CMesh &m = g.meshes[mi];

				shaderBlinnShadowMap.setupUniformsMaterial(*(m.material));

				if (m.material->texture0 != nullptr)
					m.material->texture0->bind();

				glDrawArrays(GL_TRIANGLES, m.vertex_start_id, m.vertices_count);

				if (m.material->texture0 != nullptr)
					m.material->texture0->unbind();
			}
		}


		glDisableVertexAttribArray(2);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(0);

		buffer.unbind();

		glActiveTexture(GL_TEXTURE1);
		depth_texture.setParam(GL_TEXTURE_COMPARE_MODE, GL_NONE);
		depth_texture.unbind();
		glActiveTexture(GL_TEXTURE0);

		shaderBlinnShadowMap.disable();

		CGlErrorCheck();
	}


	/**
	 * render the meshes from the object file using a shadow map texture
	 */
	void renderWithShadowAndCausticMap(	CGlLights &cGlLights,
										const GLSL::mat4 &pvm_matrix,
										const GLSL::mat4 &view_model_matrix,
										const GLSL::mat3 &view_model_normal_matrix3,
										const GLSL::mat4 &view_matrix,
										CGlTexture &depth_texture,			///< depth texture (front faces rendered from light pos)
										CGlTexture &caustic_map_texture,	///< caustic map texture (map storing the caustic light)
										CGlTexture &caustic_map_depth_texture,	///< depth texture of caustics to omit invalid caustics passing through surfaces
										const GLSL::mat4 &shadow_map_matrix
				)
	{
		shaderBlinnShadowAndCausticMap.use();

		shaderBlinnShadowAndCausticMap.setupUniformsMatrices(pvm_matrix, view_model_matrix, view_model_normal_matrix3);
		shaderBlinnShadowAndCausticMap.setupUniformsShadowMapping(shadow_map_matrix);
		shaderBlinnShadowAndCausticMap.setupUniformsLights(cGlLights, view_matrix*cGlLights.light0_world_pos4);

		glActiveTexture(GL_TEXTURE1);
		depth_texture.bind();
		depth_texture.setParam(GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
		depth_texture.setParam(GL_TEXTURE_COMPARE_FUNC, GL_GREATER);

		glActiveTexture(GL_TEXTURE2);
		caustic_map_texture.bind();

		glActiveTexture(GL_TEXTURE3);
		caustic_map_depth_texture.bind();

		glActiveTexture(GL_TEXTURE0);

		buffer.bind();

		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, buf_vertices_start_ptr);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, buf_normals_start_ptr);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, buf_texture_coords_start_ptr);
		glEnableVertexAttribArray(2);

		for (int gi = 0; gi < groups_count; gi++)
		{
			CGroup &g = groups[gi];

			for (int mi = 0; mi < g.meshes_count; mi++)
			{
				CMesh &m = g.meshes[mi];

				shaderBlinnShadowMap.setupUniformsMaterial(*(m.material));

				if (m.material->texture0 != nullptr)
					m.material->texture0->bind();

				glDrawArrays(GL_TRIANGLES, m.vertex_start_id, m.vertices_count);

				if (m.material->texture0 != nullptr)
					m.material->texture0->unbind();
			}
		}


		glDisableVertexAttribArray(2);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(0);

		buffer.unbind();

		glActiveTexture(GL_TEXTURE3);
		caustic_map_depth_texture.unbind();

		glActiveTexture(GL_TEXTURE2);
		caustic_map_texture.unbind();

		glActiveTexture(GL_TEXTURE1);
		depth_texture.setParam(GL_TEXTURE_COMPARE_MODE, GL_NONE);
		depth_texture.unbind();
		glActiveTexture(GL_TEXTURE0);

		shaderBlinnShadowAndCausticMap.disable();

		CGlErrorCheck();
	}

	void cleanup()
	{
		if (groups)
		{
			delete [] groups;
			groups = nullptr;
		}

		if (materials)
		{
			delete [] materials;
			materials = nullptr;
		}

		buffer.bind();
		buffer.data(1, nullptr);	// initialize with 1 - 0 doesn't work
		buffer.unbind();

		CGlErrorCheck();
	}

	~CGlRenderObjFile2()
	{
		cleanup();
	}
};

#endif
