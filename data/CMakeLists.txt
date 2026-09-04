if(nusyst_DOWNLOAD_DATA)
	message(STATUS "nusyst_DOWNLOAD_DATA=ON: fetching NuSystematics/nusyst_data for installation")

	# Download directly at install time so no copy is staged in build/_deps.
	install(CODE [[
		set(_dest_dir "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/data")
		set(_stage_dir "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/.nusyst_data_stage")
		set(_repo_dir "${_stage_dir}/nusyst_data-main")
		set(_repo_url "https://github.com/NuSystematics/nusyst_data.git")
		set(_repo_ref "main")
		set(_source_id "${_repo_url}#${_repo_ref}")
		set(_stamp_file "${_dest_dir}/.nusyst_data_source")

		# If already installed from the same source, skip re-downloading.
		if(EXISTS "${_stamp_file}")
			file(READ "${_stamp_file}" _installed_source)
			string(STRIP "${_installed_source}" _installed_source)
			if("${_installed_source}" STREQUAL "${_source_id}")
				message(STATUS "nusyst_data already installed in ${_dest_dir}; skipping download")
				return()
			endif()
		endif()

		message(STATUS "Downloading nusyst_data from ${_repo_url} (${_repo_ref})")
		file(MAKE_DIRECTORY "${_stage_dir}")
		find_program(_git_exe git REQUIRED)
		execute_process(
			COMMAND "${_git_exe}" clone --depth 1 --branch "${_repo_ref}" --progress "${_repo_url}" "${_repo_dir}"
			RESULT_VARIABLE _git_code
		)
		if(NOT _git_code EQUAL 0)
			message(FATAL_ERROR "Failed to clone nusyst_data repository")
		endif()

		file(MAKE_DIRECTORY "${_dest_dir}")
		execute_process(
			COMMAND "${CMAKE_COMMAND}" -E copy_directory
				"${_repo_dir}"
				"${_dest_dir}"
			RESULT_VARIABLE _copy_result
		)
		if(NOT _copy_result EQUAL 0)
			message(FATAL_ERROR "Failed to copy nusyst_data into ${_dest_dir}")
		endif()

		file(WRITE "${_stamp_file}" "${_source_id}\n")

		file(REMOVE_RECURSE "${_stage_dir}")
	]])
else()
	message(STATUS "nusyst_DOWNLOAD_DATA=OFF: installing local data directory contents")

	install(
		DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
		DESTINATION ${CMAKE_INSTALL_PREFIX}/data
		FILES_MATCHING PATTERN "*"
		PATTERN "CMakeLists.txt" EXCLUDE
	)
endif()