<template>
  <el-row class="top-nav">
    <el-menu default-active="1" class="el-menu-demo" mode="horizontal" @select="handleSelect">
      <el-col :span="22" :offset="1">
        
          <router-link to="/home"  v-if="'logo' in top">
            <el-menu-item index="1"><img class="nsp-logo" :src="'static/'+top.logo" /></el-menu-item>
          </router-link>
          <router-link to="/doc">
            <el-menu-item index="1">Documentation</el-menu-item>
          </router-link>
          <router-link to="/api">
            <el-menu-item index="1">API</el-menu-item>
          </router-link>
          <router-link to="/download">
            <el-menu-item index="1">Download</el-menu-item>
          </router-link>
          <a href="mailto:wj_hust08@hust.edu.cn?cc=yxiao@hust.edu.cn">
            <el-menu-item index="1">Contact
          </el-menu-item></a>
      </el-col>
    </el-menu>
    <div class="line"></div>
  </el-row>
</template>

<script>
  export default {
    name: 'nav',

    data() {
      return {
        top: []
      }
    },

    created() {
      this.fetch_data()
    },

    methods: {
      fetch_data() {
        var v = this
        this.$http.get('static/top.json').then((response) => {
          v.top = response.body
        }, (response) => {
        })
      },

      handleSelect(key, keyPath) {
        // console.log(key, keyPath)
      }
    }
  }
</script>

<style scoped>
  a {
    color: rgb(71, 86, 105);
    text-decoration: none;
  }

  .top-nav {
    border-bottom: 1px solid #f0f0f0;
  }

  .el-menu {
    background-color: #fff;
  }

  .nsp-logo {
    vertical-align: middle;
    height: 50px;
  }

</style>
