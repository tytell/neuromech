<?xml version='1.0' encoding='UTF-8'?>
<Project Type="Project" LVVersion="12008004">
	<Item Name="My Computer" Type="My Computer">
		<Property Name="NI.SortType" Type="Int">3</Property>
		<Property Name="server.app.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.control.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.tcp.enabled" Type="Bool">false</Property>
		<Property Name="server.tcp.port" Type="Int">0</Property>
		<Property Name="server.tcp.serviceName" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.tcp.serviceName.default" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.vi.callsEnabled" Type="Bool">true</Property>
		<Property Name="server.vi.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="specify.custom.address" Type="Bool">false</Property>
		<Item Name="TDMS to HDF5.vi" Type="VI" URL="../TDMS to HDF5.vi"/>
		<Item Name="Simple OpenCreateReplace Dataset.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/Simple OpenCreateReplace Dataset.vi"/>
		<Item Name="TDMSdatatype.ctl" Type="VI" URL="../TDMSdatatype.ctl"/>
		<Item Name="Simple H5Dwrite.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/Simple H5Dwrite.vi"/>
		<Item Name="Append Element to Dataset.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/Append Element to Dataset.vi"/>
		<Item Name="Populate HDF5 Tree.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/utility.llb/Populate HDF5 Tree.vi"/>
		<Item Name="OpenCreate Dataset.vi" Type="VI" URL="../OpenCreate Dataset.vi"/>
		<Item Name="TDMS chan to HDF5 (i16).vi" Type="VI" URL="../TDMS chan to HDF5 (i16).vi"/>
		<Item Name="TDMS chan to HDF5 (dbl).vi" Type="VI" URL="../TDMS chan to HDF5 (dbl).vi"/>
		<Item Name="TDMS chan to HDF5.vi" Type="VI" URL="../TDMS chan to HDF5.vi"/>
		<Item Name="TDMS chan to HDF5 (i32).vi" Type="VI" URL="../TDMS chan to HDF5 (i32).vi"/>
		<Item Name="TDMS chan to HDF5 (u64).vi" Type="VI" URL="../TDMS chan to HDF5 (u64).vi"/>
		<Item Name="Dependencies" Type="Dependencies">
			<Item Name="vi.lib" Type="Folder">
				<Item Name="subFile Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/express/express input/FileDialogBlock.llb/subFile Dialog.vi"/>
				<Item Name="OpenCreateReplace HDF5 File.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/file.llb/OpenCreateReplace HDF5 File.vi"/>
				<Item Name="HDF5 Ref.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/library.llb/HDF5 Ref.ctl"/>
				<Item Name="Simple Error Handler.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Simple Error Handler.vi"/>
				<Item Name="OpenCreateGroup.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/OpenCreateGroup.vi"/>
				<Item Name="H5Gclose.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/H5Gclose.vi"/>
				<Item Name="OpenCreateGroup (String).vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/OpenCreateGroup (String).vi"/>
				<Item Name="OpenCreateReplace Dataset.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/OpenCreateReplace Dataset.vi"/>
				<Item Name="LVH5Pget_default.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/property.llb/LVH5Pget_default.vi"/>
				<Item Name="H5Dopen.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/H5Dopen.vi"/>
				<Item Name="H5Pset_chunk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/property.llb/H5Pset_chunk.vi"/>
				<Item Name="H5Sclose.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sclose.vi"/>
				<Item Name="H5Tclose.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/datatype.llb/H5Tclose.vi"/>
				<Item Name="Create Matching Dataspace.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/Create Matching Dataspace.vi"/>
				<Item Name="Create Matching Datatype.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/datatype.llb/Create Matching Datatype.vi"/>
				<Item Name="H5Pcopy.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/property.llb/H5Pcopy.vi"/>
				<Item Name="H5Pset_layout.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/property.llb/H5Pset_layout.vi"/>
				<Item Name="H5D_layout_t.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/property.llb/H5D_layout_t.ctl"/>
				<Item Name="H5Pclose.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/property.llb/H5Pclose.vi"/>
				<Item Name="H5Dwrite.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/H5Dwrite.vi"/>
				<Item Name="H5Dclose.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/H5Dclose.vi"/>
				<Item Name="H5Dget_space.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/H5Dget_space.vi"/>
				<Item Name="H5Sget_simple_extent_dims.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sget_simple_extent_dims.vi"/>
				<Item Name="H5Dextend.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/H5Dextend.vi"/>
				<Item Name="H5Sselect_hyperslab.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sselect_hyperslab.vi"/>
				<Item Name="H5S_seloper_t.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5S_seloper_t.ctl"/>
				<Item Name="H5Fclose.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/file.llb/H5Fclose.vi"/>
				<Item Name="H5Pset_deflate.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/property.llb/H5Pset_deflate.vi"/>
				<Item Name="Simple H5Awrite.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/attribute.llb/Simple H5Awrite.vi"/>
				<Item Name="H5Fopen.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/file.llb/H5Fopen.vi"/>
				<Item Name="H5ErrorHandler.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/H5ErrorHandler.vi"/>
				<Item Name="h5helper.dll" Type="Document" URL="/&lt;vilib&gt;/addons/_hdf5/h5helper.dll"/>
				<Item Name="HDF5 Error Cluster.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/HDF5 Error Cluster.ctl"/>
				<Item Name="Error Cluster From Error Code.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Error Cluster From Error Code.vi"/>
				<Item Name="From HDF5 Ref.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/From HDF5 Ref.vi"/>
				<Item Name="To HDF5 Ref.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/To HDF5 Ref.vi"/>
				<Item Name="File Exists__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/File Exists__ogtk.vi"/>
				<Item Name="H5Fcreate.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/file.llb/H5Fcreate.vi"/>
				<Item Name="Not an HDF5 Refnum Constant.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/identifier.llb/Not an HDF5 Refnum Constant.vi"/>
				<Item Name="OpenCreateGroup (Array).vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/OpenCreateGroup (Array).vi"/>
				<Item Name="H5Gopen.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/H5Gopen.vi"/>
				<Item Name="H5Gcreate.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/H5Gcreate.vi"/>
				<Item Name="String to 1D Array__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/String to 1D Array__ogtk.vi"/>
				<Item Name="H5Awrite.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/attribute.llb/H5Awrite.vi"/>
				<Item Name="H5ErrorHandler2.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/H5ErrorHandler2.vi"/>
				<Item Name="H5Aclose.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/attribute.llb/H5Aclose.vi"/>
				<Item Name="H5Idec_ref.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/identifier.llb/H5Idec_ref.vi"/>
				<Item Name="LVH5Tcreate_type.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/datatype.llb/LVH5Tcreate_type.vi"/>
				<Item Name="Get TDEnum from Data__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get TDEnum from Data__ogtk.vi"/>
				<Item Name="Type Descriptor Enumeration__ogtk.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Type Descriptor Enumeration__ogtk.ctl"/>
				<Item Name="Get Header from TD__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get Header from TD__ogtk.vi"/>
				<Item Name="Type Descriptor Header__ogtk.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Type Descriptor Header__ogtk.ctl"/>
				<Item Name="Type Descriptor__ogtk.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Type Descriptor__ogtk.ctl"/>
				<Item Name="Get Array Element Default Data__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get Array Element Default Data__ogtk.vi"/>
				<Item Name="Get Array Element TD__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get Array Element TD__ogtk.vi"/>
				<Item Name="Variant to Header Info__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Variant to Header Info__ogtk.vi"/>
				<Item Name="Get Element TD from Array TD__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get Element TD from Array TD__ogtk.vi"/>
				<Item Name="Get Default Data from TD__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get Default Data from TD__ogtk.vi"/>
				<Item Name="Array of VData to VCluster__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Array of VData to VCluster__ogtk.vi"/>
				<Item Name="Set Data Name__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Set Data Name__ogtk.vi"/>
				<Item Name="Get PString__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get PString__ogtk.vi"/>
				<Item Name="Get Last PString__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get Last PString__ogtk.vi"/>
				<Item Name="Split Cluster TD__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Split Cluster TD__ogtk.vi"/>
				<Item Name="Get Data Name from TD__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get Data Name from TD__ogtk.vi"/>
				<Item Name="Array Size(s)__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Array Size(s)__ogtk.vi"/>
				<Item Name="H5Screate_simple.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Screate_simple.vi"/>
				<Item Name="To U64.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/To U64.vi"/>
				<Item Name="I32 to U64.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/I32 to U64.vi"/>
				<Item Name="I32 to U64 Array.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/I32 to U64 Array.vi"/>
				<Item Name="U32 to U64.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/U32 to U64.vi"/>
				<Item Name="I32 to U64 2D-Array.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/I32 to U64 2D-Array.vi"/>
				<Item Name="H5Screate.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Screate.vi"/>
				<Item Name="H5S_class_t.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5S_class_t.ctl"/>
				<Item Name="CreateReplace Attribute.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/attribute.llb/CreateReplace Attribute.vi"/>
				<Item Name="Check Attribute Existence.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/attribute.llb/Check Attribute Existence.vi"/>
				<Item Name="H5Aopen_name.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/attribute.llb/H5Aopen_name.vi"/>
				<Item Name="H5Acreate.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/attribute.llb/H5Acreate.vi"/>
				<Item Name="H5Adelete.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/attribute.llb/H5Adelete.vi"/>
				<Item Name="H5Dcreate.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataset.llb/H5Dcreate.vi"/>
				<Item Name="H5Gunlink.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/H5Gunlink.vi"/>
				<Item Name="Check Object Existence.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/Check Object Existence.vi"/>
				<Item Name="H5Sget_simple_extent_ndims.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sget_simple_extent_ndims.vi"/>
				<Item Name="From U64.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/From U64.vi"/>
				<Item Name="U64 to U32.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/U64 to U32.vi"/>
				<Item Name="U64 to I32.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/U64 to I32.vi"/>
				<Item Name="U64 to U32 Array.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/U64 to U32 Array.vi"/>
				<Item Name="U64 to I32 Array.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/U64 to I32 Array.vi"/>
				<Item Name="U64 to I32 2D Array.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/U64 to I32 2D Array.vi"/>
				<Item Name="From I64.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/From I64.vi"/>
				<Item Name="I64 to I32.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/I64 to I32.vi"/>
				<Item Name="I64 to I32 Array.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/I64 to I32 Array.vi"/>
				<Item Name="H5S_sel_type.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5S_sel_type.ctl"/>
				<Item Name="H5Sget_select_elem_pointlist.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sget_select_elem_pointlist.vi"/>
				<Item Name="H5Sget_select_elem_npoints.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sget_select_elem_npoints.vi"/>
				<Item Name="H5Sget_select_hyper_blocklist.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sget_select_hyper_blocklist.vi"/>
				<Item Name="H5Sget_select_hyper_nblocks.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sget_select_hyper_nblocks.vi"/>
				<Item Name="Hyperslab Blocks to Pointlist.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/Hyperslab Blocks to Pointlist.vi"/>
				<Item Name="Index Array__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Index Array__ogtk.vi"/>
				<Item Name="Reshape Array to 1D VArray__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Reshape Array to 1D VArray__ogtk.vi"/>
				<Item Name="Get Data Name__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Get Data Name__ogtk.vi"/>
				<Item Name="Compute 1D Index__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Compute 1D Index__ogtk.vi"/>
				<Item Name="H5Sget_select_type.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sget_select_type.vi"/>
				<Item Name="Array of VData to VArray__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Array of VData to VArray__ogtk.vi"/>
				<Item Name="Parse Dataset Path.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/Parse Dataset Path.vi"/>
				<Item Name="Reshape 1D Array__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Reshape 1D Array__ogtk.vi"/>
				<Item Name="Scalar to Variant Array.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/common.llb/Scalar to Variant Array.vi"/>
				<Item Name="H5Sselect_all.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/dataspace.llb/H5Sselect_all.vi"/>
				<Item Name="Not an HDF5 Refnum.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/identifier.llb/Not an HDF5 Refnum.vi"/>
				<Item Name="H5Iget_type.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/identifier.llb/H5Iget_type.vi"/>
				<Item Name="List Group Objects (recursive).vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/List Group Objects (recursive).vi"/>
				<Item Name="List Group Objects.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/List Group Objects.vi"/>
				<Item Name="H5Gget_num_objs.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/H5Gget_num_objs.vi"/>
				<Item Name="H5G_obj_t.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/H5G_obj_t.ctl"/>
				<Item Name="H5Gget_objtype_by_idx.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/H5Gget_objtype_by_idx.vi"/>
				<Item Name="H5Gget_objname_by_idx.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/group.llb/H5Gget_objname_by_idx.vi"/>
				<Item Name="Merge Errors.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Merge Errors.vi"/>
				<Item Name="compatFileDialog.vi" Type="VI" URL="/&lt;vilib&gt;/_oldvers/_oldvers.llb/compatFileDialog.vi"/>
				<Item Name="GetHelpDir.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/GetHelpDir.vi"/>
				<Item Name="BuildHelpPath.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/BuildHelpPath.vi"/>
				<Item Name="LVBoundsTypeDef.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/miscctls.llb/LVBoundsTypeDef.ctl"/>
				<Item Name="Get String Text Bounds.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Get String Text Bounds.vi"/>
				<Item Name="Get Text Rect.vi" Type="VI" URL="/&lt;vilib&gt;/picture/picture.llb/Get Text Rect.vi"/>
				<Item Name="Convert property node font to graphics font.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Convert property node font to graphics font.vi"/>
				<Item Name="Longest Line Length in Pixels.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Longest Line Length in Pixels.vi"/>
				<Item Name="Three Button Dialog CORE.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Three Button Dialog CORE.vi"/>
				<Item Name="Three Button Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Three Button Dialog.vi"/>
				<Item Name="DialogTypeEnum.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/DialogTypeEnum.ctl"/>
				<Item Name="Not Found Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Not Found Dialog.vi"/>
				<Item Name="Set Bold Text.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Set Bold Text.vi"/>
				<Item Name="Clear Errors.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Clear Errors.vi"/>
				<Item Name="eventvkey.ctl" Type="VI" URL="/&lt;vilib&gt;/event_ctls.llb/eventvkey.ctl"/>
				<Item Name="ErrWarn.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/ErrWarn.ctl"/>
				<Item Name="Details Display Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Details Display Dialog.vi"/>
				<Item Name="Search and Replace Pattern.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Search and Replace Pattern.vi"/>
				<Item Name="Find Tag.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Find Tag.vi"/>
				<Item Name="Format Message String.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Format Message String.vi"/>
				<Item Name="whitespace.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/whitespace.ctl"/>
				<Item Name="Trim Whitespace.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Trim Whitespace.vi"/>
				<Item Name="Error Code Database.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Error Code Database.vi"/>
				<Item Name="GetRTHostConnectedProp.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/GetRTHostConnectedProp.vi"/>
				<Item Name="Set String Value.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Set String Value.vi"/>
				<Item Name="TagReturnType.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/TagReturnType.ctl"/>
				<Item Name="Check Special Tags.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Check Special Tags.vi"/>
				<Item Name="General Error Handler CORE.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/General Error Handler CORE.vi"/>
				<Item Name="DialogType.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/DialogType.ctl"/>
				<Item Name="General Error Handler.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/General Error Handler.vi"/>
				<Item Name="ex_CorrectErrorChain.vi" Type="VI" URL="/&lt;vilib&gt;/express/express shared/ex_CorrectErrorChain.vi"/>
				<Item Name="Build Error Cluster__ogtk.vi" Type="VI" URL="/&lt;vilib&gt;/addons/_hdf5/openg.llb/Build Error Cluster__ogtk.vi"/>
			</Item>
			<Item Name="user.lib" Type="Folder">
				<Item Name="Conditional Auto-Indexing Tunnel__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/array/array.llb/Conditional Auto-Indexing Tunnel__ogtk.vi"/>
				<Item Name="Conditional Auto-Indexing Tunnel (I32)__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/array/array.llb/Conditional Auto-Indexing Tunnel (I32)__ogtk.vi"/>
				<Item Name="Conditional Auto-Indexing Tunnel (String)__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/array/array.llb/Conditional Auto-Indexing Tunnel (String)__ogtk.vi"/>
				<Item Name="Build Error Cluster__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/error/error.llb/Build Error Cluster__ogtk.vi"/>
			</Item>
			<Item Name="WriteConvert Multiple Attributes.vi" Type="VI" URL="../WriteConvert Multiple Attributes.vi"/>
			<Item Name="TDMS chan to HDF5 (str).vi" Type="VI" URL="../TDMS chan to HDF5 (str).vi"/>
			<Item Name="hdf5dll.dll" Type="Document" URL="hdf5dll.dll">
				<Property Name="NI.PreserveRelativePath" Type="Bool">true</Property>
			</Item>
		</Item>
		<Item Name="Build Specifications" Type="Build"/>
	</Item>
</Project>
