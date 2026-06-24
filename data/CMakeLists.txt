if(nusyst_DOWNLOAD_DATA)
	message(STATUS "nusyst_DOWNLOAD_DATA=ON: fetching NuSystematics/nusyst_data for installation")

	# Download and unpack directly at install time so no copy is staged in
	# build/_deps. We use GitHub's source tarball endpoint for the main branch.
	install(CODE [[
		set(_dest_dir "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/data")
		set(_stage_dir "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/.nusyst_data_stage")
		set(_archive "${_stage_dir}/nusyst_data-main.tar.gz")
		set(_source_url "https://codeload.github.com/NuSystematics/nusyst_data/tar.gz/refs/heads/main")
		set(_stamp_file "${_dest_dir}/.nusyst_data_source")

		# If already installed from the same source, skip re-downloading.
		if(EXISTS "${_stamp_file}")
			file(READ "${_stamp_file}" _installed_source)
			string(STRIP "${_installed_source}" _installed_source)
			if("${_installed_source}" STREQUAL "${_source_url}")
				message(STATUS "nusyst_data already installed in ${_dest_dir}; skipping download")
				return()
			endif()
		endif()

		message(STATUS "Downloading nusyst_data from ${_source_url}")
		file(MAKE_DIRECTORY "${_stage_dir}")
		file(DOWNLOAD "${_source_url}" "${_archive}" SHOW_PROGRESS STATUS _dl_status)
		list(GET _dl_status 0 _dl_code)
		list(GET _dl_status 1 _dl_msg)
		if(NOT _dl_code EQUAL 0)
			message(FATAL_ERROR "Failed to download nusyst_data: ${_dl_msg}")
		endif()

		execute_process(
			COMMAND "${CMAKE_COMMAND}" -E tar xzf "${_archive}"
			WORKING_DIRECTORY "${_stage_dir}"
			RESULT_VARIABLE _tar_result
		)
		if(NOT _tar_result EQUAL 0)
			message(FATAL_ERROR "Failed to extract nusyst_data archive")
		endif()

		file(MAKE_DIRECTORY "${_dest_dir}")
		execute_process(
			COMMAND "${CMAKE_COMMAND}" -E copy_directory
				"${_stage_dir}/nusyst_data-main"
				"${_dest_dir}"
			RESULT_VARIABLE _copy_result
		)
		if(NOT _copy_result EQUAL 0)
			message(FATAL_ERROR "Failed to copy nusyst_data into ${_dest_dir}")
		endif()

		file(WRITE "${_stamp_file}" "${_source_url}\n")

		file(REMOVE_RECURSE "${_stage_dir}")
	]])
else()
	message(STATUS "nusyst_DOWNLOAD_DATA=OFF: installing local data directory contents")

	install(
		DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
		DESTINATION ${CMAKE_INSTALL_PREFIX}/data
	)
endif()
