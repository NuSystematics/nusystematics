if(nusyst_DOWNLOAD_DATA)
	message(STATUS "nusyst_DOWNLOAD_DATA=ON: fetching NuSystematics/nusyst_data release assets for installation")

	# Download directly at install time so no copy is staged in build/_deps.
	# nusyst_data no longer stores its *.root files in git (GitHub's 100MB
	# per-file push limit, and git's inability to drop a superseded version
	# of a binary blob, made that unworkable at this data set's size) -- the
	# data lives on a GitHub Release instead, and manifest.json (checked
	# into nusyst_data's main branch) describes exactly which release, which
	# asset files, their checksums, and where each one extracts to. So this
	# is now a set of plain HTTPS downloads with hash verification, not a
	# git clone.
	install(CODE [[
		set(_dest_dir "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/data")
		set(_stage_dir "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/.nusyst_data_stage")
		set(_repo "NuSystematics/nusyst_data")
		set(_manifest_url "https://raw.githubusercontent.com/${_repo}/main/manifest.json")
		set(_stamp_file "${_dest_dir}/.nusyst_data_source")

		message(STATUS "Fetching data manifest from ${_manifest_url}")
		file(MAKE_DIRECTORY "${_stage_dir}")
		file(DOWNLOAD "${_manifest_url}" "${_stage_dir}/manifest.json"
			STATUS _dl_status
		)
		list(GET _dl_status 0 _dl_code)
		if(NOT _dl_code EQUAL 0)
			list(GET _dl_status 1 _dl_msg)
			message(FATAL_ERROR "Failed to download nusyst_data manifest.json: ${_dl_msg}")
		endif()

		file(READ "${_stage_dir}/manifest.json" _manifest_json)
		string(JSON _release_tag GET "${_manifest_json}" releaseTag)
		set(_source_id "${_repo}@${_release_tag}")

		# If already installed from the same release, skip re-downloading.
		if(EXISTS "${_stamp_file}")
			file(READ "${_stamp_file}" _installed_source)
			string(STRIP "${_installed_source}" _installed_source)
			if("${_installed_source}" STREQUAL "${_source_id}")
				message(STATUS "nusyst_data ${_release_tag} already installed in ${_dest_dir}; skipping download")
				file(REMOVE_RECURSE "${_stage_dir}")
				return()
			endif()
		endif()

		message(STATUS "Installing nusyst_data release ${_release_tag}")
		file(MAKE_DIRECTORY "${_dest_dir}")

		string(JSON _n_assets LENGTH "${_manifest_json}" assets)
		if(_n_assets GREATER 0)
			math(EXPR _last_idx "${_n_assets} - 1")
			foreach(_i RANGE ${_last_idx})
				string(JSON _asset_file   GET "${_manifest_json}" assets ${_i} file)
				string(JSON _asset_dest   GET "${_manifest_json}" assets ${_i} extract_to)
				string(JSON _asset_sha256 GET "${_manifest_json}" assets ${_i} sha256)

				set(_asset_url "https://github.com/${_repo}/releases/download/${_release_tag}/${_asset_file}")
				set(_asset_path "${_stage_dir}/${_asset_file}")

				message(STATUS "  downloading ${_asset_file}")
				file(DOWNLOAD "${_asset_url}" "${_asset_path}"
					EXPECTED_HASH SHA256=${_asset_sha256}
					STATUS _dl_status
				)
				list(GET _dl_status 0 _dl_code)
				if(NOT _dl_code EQUAL 0)
					list(GET _dl_status 1 _dl_msg)
					message(FATAL_ERROR "Failed to download ${_asset_url}: ${_dl_msg}")
				endif()

				set(_extract_dir "${_dest_dir}/${_asset_dest}")
				file(MAKE_DIRECTORY "${_extract_dir}")
				file(ARCHIVE_EXTRACT
					INPUT "${_asset_path}"
					DESTINATION "${_extract_dir}"
				)
			endforeach()
		endif()

		file(WRITE "${_stamp_file}" "${_source_id}\n")
		file(REMOVE_RECURSE "${_stage_dir}")
	]])
else()
	message(STATUS "nusyst_DOWNLOAD_DATA=OFF: installing local data directory contents")

	install(
		DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
		DESTINATION ${CMAKE_INSTALL_PREFIX}/data
	)
endif()
